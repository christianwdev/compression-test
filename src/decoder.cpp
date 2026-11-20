#include "common.h"
namespace gt {

struct DecTables {
  u32 lit[MAXCL][SCALE];
  u32 ll[SCALE], mlk[2][SCALE], off[4][SCALE], align[SCALE];
  u8 map[MAXCTX];
};

static void build_dec(const u16* f, int n, u32* tab) {
  u32 cum = 0;
  for (int s = 0; s < n; s++) {
    u32 fr = f[s];
    for (u32 j = 0; j < fr; j++) tab[cum + j] = (u32)s | (fr << 8) | (j << 20);
    cum += fr;
  }
}

static u32 LL_BASE[LL_SYMS], LL_NB[LL_SYMS], ML_BASE[MLC_SYMS], ML_NB[MLC_SYMS];
static u8 MLK_KIND[MLK_SYMS], MLK_MC[MLK_SYMS];
static bool g_init = false;
static void init_static() {
  if (g_init) return;
  for (u32 s = 0; s < LL_SYMS; s++) {
    if (s < 16) { LL_BASE[s] = s; LL_NB[s] = 0; }
    else { u32 n = 4 + (s - 16) / 2, b = (s - 16) & 1; LL_BASE[s] = (1u << n) | (b << (n - 1)); LL_NB[s] = n - 1; }
  }
  for (u32 s = 0; s < MLC_SYMS; s++) {
    if (s < 32) { ML_BASE[s] = s; ML_NB[s] = 0; }
    else { u32 n = 5 + (s - 32) / 2, b = (s - 32) & 1; ML_BASE[s] = (1u << n) | (b << (n - 1)); ML_NB[s] = n - 1; }
  }
  for (u32 s = 0; s < MLK_SYMS; s++) { MLK_KIND[s] = s / MLC_SYMS; MLK_MC[s] = s % MLC_SYMS; }
  g_init = true;
}

#define RANS_DEC(x, tab, sym) do { u32 e_ = (tab)[(x) & (SCALE - 1)]; \
    (x) = ((e_ >> 8) & 0xfff) * ((x) >> SCALE_BITS) + (e_ >> 20); (sym) = e_ & 0xff; \
    if ((x) < RANS_L) { (x) = ((x) << 16) | load16(rp); rp += 2; } } while (0)

template <int MODE>
static void decode_seqs(const DecTables& T, u8* base, u32 pos, u32 end, u32 nseq, u32* reps,
                        const u8* rp, const u8* bitp) {
  u32 x0 = load32(rp), x1 = load32(rp + 4), x2 = load32(rp + 8), x3 = load32(rp + 12); rp += 16;
  BitReader br; br.init(bitp);
  u32 r0 = reps[0], r1 = reps[1], r2 = reps[2];
  u32 li = 0;
  for (u32 s = 0; s <= nseq; s++) {
    u32 ll;
    if (s < nseq) {
      u32 sym; RANS_DEC(x2, T.ll, sym);
      ll = LL_BASE[sym] + br.get(LL_NB[sym]);
    } else ll = end - pos;
    for (u32 k = 0; k < ll; k++) {
      const u8* q = base + pos; u8 p1 = *(q - 1), mb = *(q - r0);
      const u32* t = T.lit[T.map[lit_ctx(MODE, pos, p1, mb)]];
      u32 sym;
      if (li & 1) RANS_DEC(x1, t, sym); else RANS_DEC(x0, t, sym);
      li++;
      base[pos++] = lit_sub(MODE) ? (u8)(sym + mb) : (u8)sym;
    }
    if (s == nseq) break;
    u32 sym; RANS_DEC(x3, T.mlk[ll != 0], sym);
    u32 kind = MLK_KIND[sym], mc = MLK_MC[sym];
    u32 ml = MIN_MATCH + ML_BASE[mc] + br.get(ML_NB[mc]);
    u32 off;
    if (kind == 3) {
      u32 c = ml - MIN_MATCH; if (c > 3) c = 3;
      u32 slot; RANS_DEC(x2, T.off[c], slot);
      u32 v;
      if (slot < 4) v = slot;
      else {
        u32 nb = (slot >> 1) - 1; v = (2 | (slot & 1)) << nb;
        if (nb >= 4) { v += br.get(nb - 4) << 4; u32 a; RANS_DEC(x3, T.align, a); v += a; }
        else v += br.get(nb);
      }
      off = v + 1; r2 = r1; r1 = r0; r0 = off;
    } else if (kind == 0) off = r0;
    else if (kind == 1) { off = r1; r1 = r0; r0 = off; }
    else { off = r2; r2 = r1; r1 = r0; r0 = off; }
    u8* d = base + pos; const u8* m = d - off;
    for (u32 k = 0; k < ml; k++) d[k] = m[k];
    pos += ml;
  }
  reps[0] = r0; reps[1] = r1; reps[2] = r2;
}

bool decompress(const u8* in, size_t n, std::vector<u8>& out) {
  init_static();
  if (n < 5 || memcmp(in, "GT0\1", 4)) return false;
  const u8* p = in + 4;
  u64 total = get_varint(p);
  std::vector<u8> buf(total + PAD + 64, 0);
  u8* base = buf.data() + PAD;
  DecTables* T = new DecTables;
  u32 reps[3] = {1, 2, 3};
  u64 pos = 0;
  while (pos < total) {
    u32 raw = (u32)get_varint(p);
    u8 type = *p++;
    if (type == 0) { memcpy(base + pos, p, raw); p += raw; pos += raw; continue; }
    u32 nseq = (u32)get_varint(p);
    u32 hdr = (u32)get_varint(p), rb = (u32)get_varint(p), bb = (u32)get_varint(p);
    BitReader br; br.init(p);
    int mode = br.get(4), ncl = br.get(7);
    int nctx = LIT_NCTX[mode];
    int mb = 0; while ((1 << mb) < ncl) mb++;
    for (int i = 0; i < nctx; i++) T->map[i] = (u8)br.get(mb);
    u16 f[256];
    for (int c = 0; c < ncl; c++) { read_table(br, f, 256); build_dec(f, 256, T->lit[c]); }
    read_table(br, f, LL_SYMS); build_dec(f, LL_SYMS, T->ll);
    for (int c = 0; c < 2; c++) { read_table(br, f, MLK_SYMS); build_dec(f, MLK_SYMS, T->mlk[c]); }
    for (int c = 0; c < 4; c++) { read_table(br, f, OFF_SYMS); build_dec(f, OFF_SYMS, T->off[c]); }
    read_table(br, f, ALIGN_SYMS); build_dec(f, ALIGN_SYMS, T->align);
    const u8* rp = p + hdr; const u8* bitp = rp + rb;
    u32 end = (u32)(pos + raw);
    switch (mode) {
      case 0: decode_seqs<0>(*T, base, (u32)pos, end, nseq, reps, rp, bitp); break;
      case 1: decode_seqs<1>(*T, base, (u32)pos, end, nseq, reps, rp, bitp); break;
      case 2: decode_seqs<2>(*T, base, (u32)pos, end, nseq, reps, rp, bitp); break;
      case 3: decode_seqs<3>(*T, base, (u32)pos, end, nseq, reps, rp, bitp); break;
      case 4: decode_seqs<4>(*T, base, (u32)pos, end, nseq, reps, rp, bitp); break;
      case 5: decode_seqs<5>(*T, base, (u32)pos, end, nseq, reps, rp, bitp); break;
      default: decode_seqs<6>(*T, base, (u32)pos, end, nseq, reps, rp, bitp); break;
    }
    p = bitp + bb; pos += raw;
  }
  delete T;
  out.assign(base, base + total);
  return true;
}

}  // namespace gt
