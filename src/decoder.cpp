#include "common.h"
#ifdef __linux__
#include <sys/mman.h>
#endif
namespace gt {

struct DecTables {
  u32 lit[MAXCL][1 << LIT_BITS];
  u32 ll[LL_MAXCL][1 << LL_BITS], mlk[MLK_MAXCL][1 << MLK_BITS], off[OFF_MAXCL][1 << OFF_BITS], align[AL_MAXCL][1 << AL_BITS];
  u8 litmap[NT * MAXCTX], llmap[NT * LL_NCTX], mlkmap[NT * MLK_NCTX], offmap[NT * OFF_NCTX], almap[NT];
  u8 ltype[MAXCH], ctype[MAXCH], dtype[MAXCH];
  u8 lmode[NT]; u32 lbase[NT];
  int lgchunk; bool rawalign;
};

// Returns false if the frequencies do not sum to the table size (corrupt input).
static bool build_dec(const u16* f, int n, u32* tab, int bits) {
  u32 cum = 0; const u32 size = 1u << bits;
  for (int s = 0; s < n; s++) {
    u32 fr = f[s];
    if (cum + fr > size) return false;
    u32 e = (u32)s | (fr << 8);
    for (u32 j = 0; j < fr; j++) tab[cum + j] = e | (j << 20);
    cum += fr;
  }
  return cum == size;
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

#define RANS_DEC(x, tab, sym, BITS) do { u32 e_ = (tab)[(x) & ((1u << (BITS)) - 1)]; \
    (x) = ((e_ >> 8) & 0xfff) * ((x) >> (BITS)) + (e_ >> 20); (sym) = e_ & 0xff; \
    u32 w_ = load16(rp); u32 m_ = 0u - (u32)((x) < RANS_L); (x) = ((((x) << 16) | w_) & m_) | ((x) & ~m_); rp += m_ & 2; } while (0)

// smallest multiple of the offset that is >= 16
static const u8 SMALL_OFF_PERIOD[16] = {0, 16, 16, 18, 16, 20, 18, 21, 16, 18, 20, 22, 24, 26, 28, 30};
// __builtin_memcpy: plain memcpy here becomes a checked library call under _FORTIFY_SOURCE
static inline void copy16(u8* d, const u8* s) { u64 a = load64(s), b = load64(s + 8); __builtin_memcpy(d, &a, 8); __builtin_memcpy(d + 8, &b, 8); }

static void decode_seqs(const DecTables& T, u8* base, const u32 bs, const u32 end, u32 nseq, u32* reps,
                        const u8* rp, const u8* bitp) {
  u32 xa = load32(rp), xb = load32(rp + 4), x2 = load32(rp + 8), x3 = load32(rp + 12), x4 = load32(rp + 16); rp += 4 * NUM_STREAMS;
  BitReader br; br.init(bitp);
  u32 r0 = reps[0], r1 = reps[1], r2 = reps[2], r3 = reps[3];
  u32 pos = bs, chunk_end = bs;
  const int lg = T.lgchunk;
  const u8* lmap = T.litmap; const u8* mmap = T.mlkmap; const u8* omap = T.offmap;
  const u8* llm = T.llmap; const u32* alt = T.align[0];
  const LitMode* lm = &lit_mode(0); u32 submask = 0;
  const bool rawalign = T.rawalign;
#define UPDATE_CHUNK() do { if (pos < end) { u32 ch = (pos - bs) >> lg; chunk_end = bs + ((ch + 1) << lg); \
    u32 lt_ = T.ltype[ch]; lm = &lit_mode(T.lmode[lt_]); submask = lm->sub ? 0xff : 0; lmap = T.litmap + T.lbase[lt_]; \
    llm = T.llmap + T.ctype[ch] * LL_NCTX; \
    mmap = T.mlkmap + T.ctype[ch] * MLK_NCTX; omap = T.offmap + T.dtype[ch] * OFF_NCTX; alt = T.align[T.almap[T.dtype[ch]]]; } } while (0)
#define LIT1(X) do { u32 mb_ = *(q - r0); \
    const u32* t_ = T.lit[lmap[lm->A[p1] + lm->B[mb_] + lm->C[p2] + lm->D[(size_t)(q - base) & 3]]]; \
    u32 sym_; RANS_DEC(X, t_, sym_, LIT_BITS); p2 = p1; p1 = (u8)(sym_ + (mb_ & submask)); *q++ = (u8)p1; } while (0)
  for (u32 s = 0; s <= nseq; s++) {
    u32 ll;
    if (pos >= chunk_end) UPDATE_CHUNK();
    br.refill();
    if (s < nseq) {
      u32 sym; RANS_DEC(x2, T.ll[llm[pos & 3]], sym, LL_BITS);
      u32 nb = LL_NB[sym]; ll = LL_BASE[sym] + (u32)(br.acc & ((1ull << nb) - 1)); br.acc >>= nb; br.n -= nb;
    } else ll = end - pos;
    const u32 ll0 = ll;
    if (ll) {
      u32 p1 = base[(size_t)pos - 1], p2 = base[(size_t)pos - 2];
      do {
        u32 n = std::min(ll, chunk_end - pos);
        u8* q = base + pos; u8* qe = q + n;
        while (q + 2 <= qe) { LIT1(xa); LIT1(xb); }
        if (q < qe) { LIT1(xa); u32 t = xa; xa = xb; xb = t; }
        pos += n; ll -= n;
        if (pos >= chunk_end) UPDATE_CHUNK();
      } while (ll);
    }
    if (s == nseq) break;
    u32 sym; RANS_DEC(x3, T.mlk[mmap[(pos & 3) * MLK_LLCTX + mlk_ctx(ll0)]], sym, MLK_BITS);
    u32 kind = MLK_KIND[sym], mc = MLK_MC[sym];
    u32 mnb = ML_NB[mc];
    u32 ml = MIN_MATCH + ML_BASE[mc] + (u32)(br.acc & ((1ull << mnb) - 1)); br.acc >>= mnb; br.n -= mnb;
    u32 off;
    if (kind == 3) {
      u32 c = ml - MIN_MATCH; if (c > 3) c = 3;
      u32 slot; RANS_DEC(x4, T.off[omap[(pos & 3) * OFF_LENCTX + c]], slot, OFF_BITS);
      if (slot < OFF_SLOTS) {
        u32 v;
        if (slot < 4) v = slot;
        else {
          u32 nb = (slot >> 1) - 1; v = (2 | (slot & 1)) << nb;
          if (nb >= 4 && !rawalign) { v += br.get(nb - 4) << 4; u32 a; RANS_DEC(x3, alt, a, AL_BITS); v += a; }
          else v += br.get(nb);
        }
        off = v + 1;
      } else if (slot == OFF_SLOTS) off = r3;
      else { u32 idx = slot - OFF_SLOTS - 1; off = (idx < 6 ? r0 : r1) + NEAR_DELTA[idx % 6]; }
      r3 = r2; r2 = r1; r1 = r0; r0 = off;
    } else if (kind == 0) off = r0;
    else if (kind == 1) { off = r1; r1 = r0; r0 = off; }
    else { off = r2; r2 = r1; r1 = r0; r0 = off; }
    u8* d = base + pos; const u8* m = d - off; u8* de = d + ml;
    if (off >= 16) {
      copy16(d, m); copy16(d + 16, m + 16);
      if (ml > 32) { d += 32; m += 32; do { copy16(d, m); d += 16; m += 16; } while (d < de); }
    } else {
      // overlapping copy: lay down 16 bytes one at a time, then continue with a period that is a
      // multiple of the offset and at least 16 so that 16-byte copies are safe
      for (int k = 0; k < 16; k++) d[k] = m[k];
      if (ml > 16) { m = d + 16 - SMALL_OFF_PERIOD[off]; d += 16; do { copy16(d, m); d += 16; m += 16; } while (d < de); }
    }
    pos += ml;
  }
  reps[0] = r0; reps[1] = r1; reps[2] = r2; reps[3] = r3;
#undef UPDATE_CHUNK
#undef LIT1
}

// Buffered input: need(n) returns a pointer to n contiguous unread bytes (plus 8 readable bytes of
// slack for the bit reader), or nullptr at end of input.
struct In {
  FILE* f; std::vector<u8> buf; size_t pos = 0, len = 0;
  explicit In(FILE* f_) : f(f_) { buf.resize(1 << 16); }
  const u8* need(size_t n) {
    if (len - pos < n) {
      if (pos) { memmove(buf.data(), buf.data() + pos, len - pos); len -= pos; pos = 0; }
      if (buf.size() < n + 16) buf.resize(n + 16);
      while (len < n) {
        size_t got = fread(buf.data() + len, 1, buf.size() - 16 - len, f);
        if (!got) break;
        len += got;
      }
      if (len < n) return nullptr;
      memset(buf.data() + len, 0, 16);
    }
    return buf.data() + pos;
  }
  size_t avail_hint() { need(32 < len - pos ? len - pos : 32); return len - pos; }
};

// Streaming decoder: keeps one window of history plus a quarter window of headroom, so memory does
// not grow with the file.
bool decompress(FILE* fin, FILE* fout) {
  init_static();
  In in(fin);
  const u8* p = in.need(6);
  if (!p || memcmp(p, "GT0\7", 4)) return false;
  in.pos += 4;
  size_t av = in.avail_hint(); if (av < 3) return false;
  p = in.buf.data() + in.pos; const u8* p0 = p;
  u64 total = get_varint(p); int lgwin = *p++; int lgblock = *p++;
  in.pos += p - p0;
  const u64 W = ((u64)1 << lgwin) + 16, B = (u64)1 << lgblock;
  const u64 cap = std::min<u64>(total, W + std::max<u64>(W / 4, B));
  // 2 MB alignment plus a huge-page hint: the output buffer is touched once, front to back, and
  // per-4K page faults would otherwise be a visible share of decode time.
  const size_t HP = (size_t)2 << 20; size_t msize = (cap + PAD + 96 + HP - 1) & ~(HP - 1);
  u8* mem = (u8*)aligned_alloc(HP, msize); if (!mem) return false;
#ifdef __linux__
  madvise(mem, msize, MADV_HUGEPAGE);
#endif
  memset(mem, 0, PAD);
  u8* base = mem + PAD;
  DecTables* T = new DecTables;
  u32 reps[NUM_REPS] = {1, 2, 3, 4};
  u64 done = 0; u32 pos = 0;  // pos: write position inside buf
  bool ok = true;
  while (done < total) {
    in.avail_hint(); p = in.buf.data() + in.pos; p0 = p;
    if (in.len - in.pos < 2) { ok = false; break; }
    u32 raw = (u32)get_varint(p);
    u8 type = *p++;
    if (raw > B || raw > total - done) { ok = false; break; }
    if ((u64)pos + raw > cap) {  // slide the window; keep the shift a multiple of 4 (position contexts)
      u32 keep = (u32)std::min<u64>(pos, W); u32 shift = (pos - keep) & ~3u; keep = pos - shift;
      memmove(base, base + shift, keep); pos = keep;
    }
    if (type == 0) {
      in.pos += p - p0;
      const u8* d = in.need(raw); if (!d) { ok = false; break; }
      memcpy(base + pos, d, raw); in.pos += raw;
    } else {
      u32 nseq = (u32)get_varint(p);
      u32 hdr = (u32)get_varint(p), rb = (u32)get_varint(p), bb = (u32)get_varint(p);
      in.pos += p - p0;
      p = in.need((size_t)hdr + rb + bb); if (!p) { ok = false; break; }
      RcDec br; br.init(p); MapModel mm[8]; TabModel tml, tmll, tmmk, tmof, tmal;
      int lg = (int)br.direct(5); T->lgchunk = lg;
      int nch = (int)((raw + (1u << lg) - 1) >> lg);
      if (nch > MAXCH) { ok = false; break; }
      int nlt = (int)br.direct(4) + 1; read_map(br, mm[0], T->ltype, nch, nlt);
      int nct = (int)br.direct(4) + 1; read_map(br, mm[1], T->ctype, nch, nct);
      int ndt = (int)br.direct(4) + 1; read_map(br, mm[2], T->dtype, nch, ndt);
      u32 nmap = 0;
      for (int t = 0; t < nlt; t++) {
        int mode = (int)br.direct(4); if (mode >= NUM_LITMODES) { ok = false; break; }
        T->lmode[t] = (u8)mode; T->lbase[t] = nmap; nmap += lit_mode(mode).nctx;
      }
      if (!ok) break;
      int ncl = (int)br.direct(6) + 1;
      read_map(br, mm[3], T->litmap, (int)nmap, ncl);
      u16 f[256];
      for (int c = 0; c < ncl; c++) { if (!read_table(br, tml, f, 256, LIT_BITS) || !build_dec(f, 256, T->lit[c], LIT_BITS)) ok = false; }
      int nll = (int)br.direct(3) + 1; read_map(br, mm[4], T->llmap, nct * LL_NCTX, nll);
      for (int c = 0; c < nll; c++) { if (!read_table(br, tmll, f, LL_SYMS, LL_BITS) || !build_dec(f, LL_SYMS, T->ll[c], LL_BITS)) ok = false; }
      int nmk = (int)br.direct(4) + 1; read_map(br, mm[5], T->mlkmap, nct * MLK_NCTX, nmk);
      for (int c = 0; c < nmk; c++) { if (!read_table(br, tmmk, f, MLK_SYMS, MLK_BITS) || !build_dec(f, MLK_SYMS, T->mlk[c], MLK_BITS)) ok = false; }
      int nof = (int)br.direct(4) + 1; read_map(br, mm[6], T->offmap, ndt * OFF_NCTX, nof);
      for (int c = 0; c < nof; c++) { if (!read_table(br, tmof, f, OFF_SYMS, OFF_BITS) || !build_dec(f, OFF_SYMS, T->off[c], OFF_BITS)) ok = false; }
      int nal = (int)br.direct(3); T->rawalign = nal == 0;  // 0: low offset bits are stored raw
      if (nal > AL_MAXCL) { ok = false; break; }
      if (nal) read_map(br, mm[7], T->almap, ndt, nal); else memset(T->almap, 0, sizeof T->almap);
      for (int c = 0; c < nal; c++) { if (!read_table(br, tmal, f, ALIGN_SYMS, AL_BITS) || !build_dec(f, ALIGN_SYMS, T->align[c], AL_BITS)) ok = false; }
      if (!ok) break;
      const u8* rp = p + hdr; const u8* bitp = rp + rb;
      u32 end = pos + raw;
      decode_seqs(*T, base, pos, end, nseq, reps, rp, bitp);
      in.pos += (size_t)hdr + rb + bb;
    }
    if (fwrite(base + pos, 1, raw, fout) != raw) { ok = false; break; }
    pos += raw; done += raw;
  }
  delete T; free(mem);
  return ok;
}

}  // namespace gt
