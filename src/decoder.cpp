#include "common.h"
namespace gt {

struct DecTables {
  u32 lit[MAXCL][SCALE];
  u32 ll[LL_MAXCL][SCALE], mlk[MLK_MAXCL][SCALE], off[OFF_MAXCL][SCALE], align[AL_MAXCL][SCALE];
  u8 litmap[NT * MAXCTX], llmap[NT], mlkmap[NT * MLK_NCTX], offmap[NT * 4], almap[NT];
  u8 ltype[MAXCH], ctype[MAXCH], dtype[MAXCH];
  int lgchunk;
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
static void decode_seqs(const DecTables& T, u8* base, const u32 bs, const u32 end, u32 nseq, u32* reps,
                        const u8* rp, const u8* bitp) {
  u32 x0 = load32(rp), x1 = load32(rp + 4), x2 = load32(rp + 8), x3 = load32(rp + 12); rp += 16;
  BitReader br; br.init(bitp);
  u32 r0 = reps[0], r1 = reps[1], r2 = reps[2], r3 = reps[3];
  u32 li = 0, pos = bs, chunk_end = bs;
  const int lg = T.lgchunk; const int nctx = LIT_NCTX[MODE];
  const u8* lmap = T.litmap; const u8* mmap = T.mlkmap; const u8* omap = T.offmap;
  const u32* llt = T.ll[0]; const u32* alt = T.align[0];
#define UPDATE_CHUNK() do { if (pos < end) { u32 ch = (pos - bs) >> lg; chunk_end = bs + ((ch + 1) << lg); \
    lmap = T.litmap + (size_t)T.ltype[ch] * nctx; llt = T.ll[T.llmap[T.ctype[ch]]]; \
    mmap = T.mlkmap + T.ctype[ch] * MLK_NCTX; omap = T.offmap + T.dtype[ch] * 4; alt = T.align[T.almap[T.dtype[ch]]]; } } while (0)
  for (u32 s = 0; s <= nseq; s++) {
    u32 ll;
    if (pos >= chunk_end) UPDATE_CHUNK();
    if (s < nseq) {
      u32 sym; RANS_DEC(x2, llt, sym);
      ll = LL_BASE[sym] + br.get(LL_NB[sym]);
    } else ll = end - pos;
    const u32 ll0 = ll;
    while (ll) {
      u32 n = std::min(ll, chunk_end - pos);
      for (u32 k = 0; k < n; k++) {
        const u8* q = base + pos; u8 p1 = *(q - 1), mb = *(q - r0);
        const u32* t = T.lit[lmap[lit_ctx(MODE, pos, p1, *(q - 2), mb)]];
        u32 sym;
        if (li & 1) RANS_DEC(x1, t, sym); else RANS_DEC(x0, t, sym);
        li++;
        base[pos++] = lit_sub(MODE) ? (u8)(sym + mb) : (u8)sym;
      }
      ll -= n;
      if (pos >= chunk_end) UPDATE_CHUNK();
    }
    if (s == nseq) break;
    u32 sym; RANS_DEC(x3, T.mlk[mmap[mlk_ctx(ll0)]], sym);
    u32 kind = MLK_KIND[sym], mc = MLK_MC[sym];
    u32 ml = MIN_MATCH + ML_BASE[mc] + br.get(ML_NB[mc]);
    u32 off;
    if (kind == 3) {
      u32 c = ml - MIN_MATCH; if (c > 3) c = 3;
      u32 slot; RANS_DEC(x2, T.off[omap[c]], slot);
      if (slot < OFF_SLOTS) {
        u32 v;
        if (slot < 4) v = slot;
        else {
          u32 nb = (slot >> 1) - 1; v = (2 | (slot & 1)) << nb;
          if (nb >= 4) { v += br.get(nb - 4) << 4; u32 a; RANS_DEC(x3, alt, a); v += a; }
          else v += br.get(nb);
        }
        off = v + 1;
      } else if (slot == OFF_SLOTS) off = r3;
      else { u32 idx = slot - OFF_SLOTS - 1; off = (idx < 6 ? r0 : r1) + NEAR_DELTA[idx % 6]; }
      r3 = r2; r2 = r1; r1 = r0; r0 = off;
    } else if (kind == 0) off = r0;
    else if (kind == 1) { off = r1; r1 = r0; r0 = off; }
    else { off = r2; r2 = r1; r1 = r0; r0 = off; }
    u8* d = base + pos; const u8* m = d - off;
    for (u32 k = 0; k < ml; k++) d[k] = m[k];
    pos += ml;
  }
  reps[0] = r0; reps[1] = r1; reps[2] = r2; reps[3] = r3;
#undef UPDATE_CHUNK
}

bool decompress(const u8* in, size_t n, std::vector<u8>& out) {
  init_static();
  if (n < 5 || memcmp(in, "GT0\2", 4)) return false;
  const u8* p = in + 4;
  u64 total = get_varint(p);
  std::vector<u8> buf(total + PAD + 64, 0);
  u8* base = buf.data() + PAD;
  DecTables* T = new DecTables;
  u32 reps[NUM_REPS] = {1, 2, 3, 4};
  u64 pos = 0;
  while (pos < total) {
    u32 raw = (u32)get_varint(p);
    u8 type = *p++;
    if (type == 0) { memcpy(base + pos, p, raw); p += raw; pos += raw; continue; }
    u32 nseq = (u32)get_varint(p);
    u32 hdr = (u32)get_varint(p), rb = (u32)get_varint(p), bb = (u32)get_varint(p);
    BitReader br; br.init(p);
    int lg = br.get(5); T->lgchunk = lg;
    int nch = (int)((raw + (1u << lg) - 1) >> lg);
    int nlt = br.get(4) + 1; read_map(br, T->ltype, nch, nlt);
    int nct = br.get(4) + 1; read_map(br, T->ctype, nch, nct);
    int ndt = br.get(4) + 1; read_map(br, T->dtype, nch, ndt);
    int mode = br.get(4), ncl = br.get(6) + 1;
    int nctx = LIT_NCTX[mode];
    read_map(br, T->litmap, nlt * nctx, ncl);
    u16 f[256];
    for (int c = 0; c < ncl; c++) { read_table(br, f, 256); build_dec(f, 256, T->lit[c]); }
    int nll = br.get(3) + 1; read_map(br, T->llmap, nct, nll);
    for (int c = 0; c < nll; c++) { read_table(br, f, LL_SYMS); build_dec(f, LL_SYMS, T->ll[c]); }
    int nmk = br.get(4) + 1; read_map(br, T->mlkmap, nct * MLK_NCTX, nmk);
    for (int c = 0; c < nmk; c++) { read_table(br, f, MLK_SYMS); build_dec(f, MLK_SYMS, T->mlk[c]); }
    int nof = br.get(4) + 1; read_map(br, T->offmap, ndt * 4, nof);
    for (int c = 0; c < nof; c++) { read_table(br, f, OFF_SYMS); build_dec(f, OFF_SYMS, T->off[c]); }
    int nal = br.get(2) + 1; read_map(br, T->almap, ndt, nal);
    for (int c = 0; c < nal; c++) { read_table(br, f, ALIGN_SYMS); build_dec(f, ALIGN_SYMS, T->align[c]); }
    const u8* rp = p + hdr; const u8* bitp = rp + rb;
    u32 end = (u32)(pos + raw);
    switch (mode) {
      case 0: decode_seqs<0>(*T, base, (u32)pos, end, nseq, reps, rp, bitp); break;
      case 1: decode_seqs<1>(*T, base, (u32)pos, end, nseq, reps, rp, bitp); break;
      case 2: decode_seqs<2>(*T, base, (u32)pos, end, nseq, reps, rp, bitp); break;
      case 3: decode_seqs<3>(*T, base, (u32)pos, end, nseq, reps, rp, bitp); break;
      case 4: decode_seqs<4>(*T, base, (u32)pos, end, nseq, reps, rp, bitp); break;
      case 5: decode_seqs<5>(*T, base, (u32)pos, end, nseq, reps, rp, bitp); break;
      case 6: decode_seqs<6>(*T, base, (u32)pos, end, nseq, reps, rp, bitp); break;
      default: decode_seqs<7>(*T, base, (u32)pos, end, nseq, reps, rp, bitp); break;
    }
    p = bitp + bb; pos += raw;
  }
  delete T;
  out.assign(base, base + total);
  return true;
}

}  // namespace gt
