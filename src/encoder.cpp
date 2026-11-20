#include "common.h"
#include <chrono>
namespace gt {
static double g_time[6];
static inline double now() { return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count(); }

static const u32 NIL = 0xffffffffu;
static const float INF = 1e30f;

struct Match { u32 len, off; };
// Matches for one block, stored compactly: per-position count, then (offset, length) pairs.
struct MatchStore { std::vector<u8> cnt; std::vector<u32> off; std::vector<u16> len; };
struct Seq { u32 ll, ml, off, kind; };
struct Hist { u32 c[256]; u32 tot; };

static inline u32 match_len(const u8* p, const u8* q, u32 maxl) {
  u32 l = 0;
  while (l + 8 <= maxl) {
    u64 d = load64(p + l) ^ load64(q + l);
    if (d) return l + (__builtin_ctzll(d) >> 3);
    l += 8;
  }
  while (l < maxl && p[l] == q[l]) l++;
  return l;
}

// ---------------- binary-tree match finder (LZMA bt4 style) ----------------
struct MF {
  const u8* b; u32 n, W, wmask, nice; int depth, hb;
  std::vector<u32> head, son, h2, h3;
  // b is the (virtual) address of file position 0; n is the end of the data currently loaded.
  MF(u64 total, const Params& P) {
    b = nullptr; n = 0; nice = P.nice; depth = P.depth; hb = P.hashbits;
    u32 w = 1u << P.lgwin; while (w / 2 >= total && w > 65536) w /= 2;
    W = w; wmask = W - 1;
    head.assign((size_t)1 << hb, NIL); son.resize(2 * (size_t)std::min<u64>(W, total + 1));
    h2.assign(65536, NIL); h3.assign(65536, NIL);
  }
  u32 find(u32 pos, Match* out, bool want) {
    u32 avail = n - pos; if (avail < 4) return 0;
    const u8* cur = b + pos; u32 maxlen = std::min(nice, avail); u32 cnt = 0, best = 1;
    u32 v = load32(cur);
    u32 k2 = v & 0xffff, k3 = ((v & 0xffffff) * 506832829u) >> 16, k4 = (v * 2654435761u) >> (32 - hb);
    u32 c2 = h2[k2], c3 = h3[k3], c = head[k4];
    h2[k2] = h3[k3] = head[k4] = pos;
    if (want) {
      if (c2 != NIL && pos - c2 < W) { u32 l = match_len(cur, b + c2, maxlen); out[cnt++] = {l, pos - c2}; best = l; }
      if (c3 != NIL && c3 != c2 && pos - c3 < W) {
        u32 l = match_len(cur, b + c3, maxlen);
        if (l > best) { out[cnt++] = {l, pos - c3}; best = l; }
      }
    }
    u32* ptr0 = &son[2 * (size_t)(pos & wmask) + 1]; u32* ptr1 = &son[2 * (size_t)(pos & wmask)];
    u32 len0 = 0, len1 = 0; int d = depth;
    while (true) {
      if (c == NIL || pos - c >= W || d-- == 0) { *ptr0 = *ptr1 = NIL; break; }
      u32* pair = &son[2 * (size_t)(c & wmask)];
      const u8* pb = b + c; u32 len = std::min(len0, len1);
      if (pb[len] == cur[len]) {
        len += 1; len += match_len(cur + len, pb + len, maxlen - len);
        if (want && len > best) { best = len; out[cnt++] = {len, pos - c}; }
        if (len == maxlen) { *ptr1 = pair[0]; *ptr0 = pair[1]; break; }
      }
      if (pb[len] < cur[len]) { *ptr1 = c; ptr1 = pair + 1; c = *ptr1; len1 = len; }
      else { *ptr0 = c; ptr0 = pair; c = *ptr0; len0 = len; }
    }
    if (want && cnt && out[cnt - 1].len == maxlen && maxlen == nice && avail > nice)
      out[cnt - 1].len = match_len(cur, cur - out[cnt - 1].off, std::min<u32>(avail, 65535));
    return cnt;
  }
};

// ---------------- models / prices ----------------
// A block is divided into fixed-size chunks. Each chunk has a literal type, a command type and a
// distance type; (type, context) pairs map to shared entropy tables.
struct Models {
  int lgchunk, nch, nlt, nct, ndt, period;
  u8 ltype[MAXCH], ctype[MAXCH], dtype[MAXCH];
  int ncl, nll, nmk, nof, nal;
  u8 lmode[NT]; u32 lbase[NT];  // per literal type: context mode and offset into litmap
  u32 litn[MAXCL], lln[LL_MAXCL], mlkn[MLK_MAXCL], offn[OFF_MAXCL], aln[AL_MAXCL];  // symbols coded per table
  u8 litmap[NT * MAXCTX], llmap[NT * LL_MAXCTX], mlkmap[NT * MLK_MAXCTX], offmap[NT * OFF_MAXCTX], almap[NT];
  u16 lit[MAXCL][256], ll[LL_MAXCL][LL_SYMS], mlk[MLK_MAXCL][MLK_SYMS], off[OFF_MAXCL][OFF_SYMS], align[AL_MAXCL][ALIGN_SYMS];
};
struct Prices {
  int lgchunk, nch, period;
  u8 lmode[NT]; u32 lbase[NT];
  u8 ltype[MAXCH], ctype[MAXCH], dtype[MAXCH];
  u8 litmap[NT * MAXCTX], llmap[NT * LL_MAXCTX], mlkmap[NT * MLK_MAXCTX], offmap[NT * OFF_MAXCTX], almap[NT];
  float lit[MAXCL][256], ll[LL_MAXCL][LL_SYMS], mlk[MLK_MAXCL][MLK_SYMS], off[OFF_MAXCL][OFF_SYMS], align[AL_MAXCL][ALIGN_SYMS];
};
static float FPRICE[(1 << MAX_BITS) + 1];  // price of frequency f out of 1 << MAX_BITS
static double XLOGX[65536];
static bool g_einit = false;
static void einit() {
  if (g_einit) return;
  FPRICE[0] = 13.0f; for (int i = 1; i <= (1 << MAX_BITS); i++) FPRICE[i] = (float)log2((double)(1 << MAX_BITS) / i);
  XLOGX[0] = 0; for (int i = 1; i < 65536; i++) XLOGX[i] = i * log2((double)i);
  g_einit = true;
}
static inline double xlogx(u32 x) { return x < 65536 ? XLOGX[x] : x * log2((double)x); }

static void default_prices(Prices& P, int lgchunk) {
  memset(&P, 0, sizeof P);
  P.lgchunk = lgchunk; P.nch = 1; P.period = 4;
  for (int s = 0; s < 256; s++) P.lit[0][s] = 6.5f;
  for (int s = 0; s < LL_SYMS; s++) P.ll[0][s] = 2.0f + 0.25f * s;
  for (int s = 0; s < MLK_SYMS; s++) {
    int k = s / MLC_SYMS, mc = s % MLC_SYMS;
    P.mlk[0][s] = (k == 0 ? 3.0f : k == 3 ? 4.0f : 6.0f) + 0.1f * std::min(mc, 32);
  }
  for (int s = 0; s < OFF_SYMS; s++) P.off[0][s] = s < OFF_SLOTS ? 5.5f : 8.0f;
  for (int s = 0; s < ALIGN_SYMS; s++) P.align[0][s] = 4.0f;
}
static int g_coarse = 1;
static double g_unseen_cap = 14.0;
static int g_smooth = 1;
// Prices for the parser. A symbol a table has not seen is not priced at a flat maximum: it is priced
// as an escape (probability from how many distinct symbols the table holds versus how many it has
// coded) times the symbol's share across all tables of the same alphabet. Otherwise fine-grained
// contexts look far more expensive to the parser than they are once the tables are rebuilt.
template <int NS> static void smooth_prices(const u16 (*f)[NS], const u32* n, int ncl, float (*p)[NS], int bits) {
  const int SCALE = 1 << bits; const int sh = MAX_BITS - bits;
  double q[NS]; double N = 0;
  for (int s = 0; s < NS; s++) q[s] = 0;
  for (int c = 0; c < ncl; c++) { N += (double)n[c] * SCALE; for (int s = 0; s < NS; s++) q[s] += (double)n[c] * f[c][s]; }
  for (int c = 0; c < ncl; c++) {
    int nz = 0; for (int s = 0; s < NS; s++) nz += f[c][s] != 0;
    double e = (nz + 1.0) / ((double)n[c] + nz + 1.0);
    for (int s = 0; s < NS; s++) {
      if (f[c][s]) { p[c][s] = FPRICE[(u32)f[c][s] << sh]; continue; }
      double v = 13.0;
      if (g_smooth && N > 0) { double qs = std::max(q[s] / N, 1.0 / 65536); v = -log2(e * qs); if (v > g_unseen_cap) v = g_unseen_cap; if (v < 1.0) v = 1.0; }
      p[c][s] = (float)v;
    }
  }
}
static void prices_from_models(const Models& M, Prices& P) {
  P.lgchunk = M.lgchunk; P.nch = M.nch; P.period = M.period; memcpy(P.lmode, M.lmode, sizeof P.lmode); memcpy(P.lbase, M.lbase, sizeof P.lbase);
  memcpy(P.ltype, M.ltype, sizeof P.ltype); memcpy(P.ctype, M.ctype, sizeof P.ctype); memcpy(P.dtype, M.dtype, sizeof P.dtype);
  memcpy(P.litmap, M.litmap, sizeof P.litmap); memcpy(P.llmap, M.llmap, sizeof P.llmap);
  memcpy(P.mlkmap, M.mlkmap, sizeof P.mlkmap); memcpy(P.offmap, M.offmap, sizeof P.offmap); memcpy(P.almap, M.almap, sizeof P.almap);
  smooth_prices<256>(M.lit, M.litn, M.ncl, P.lit, LIT_BITS);
  smooth_prices<LL_SYMS>(M.ll, M.lln, M.nll, P.ll, LL_BITS);
  smooth_prices<MLK_SYMS>(M.mlk, M.mlkn, M.nmk, P.mlk, MLK_BITS);
  smooth_prices<OFF_SYMS>(M.off, M.offn, M.nof, P.off, OFF_BITS);
  if (M.nal) smooth_prices<ALIGN_SYMS>(M.align, M.aln, M.nal, P.align, AL_BITS);
  else { memset(P.almap, 0, sizeof P.almap); for (int s = 0; s < ALIGN_SYMS; s++) P.align[0][s] = 4.0f; }
}
// Reuse the previous block's statistics as the starting guess for a new block: every chunk of the
// new block takes the types of the previous block's last chunk.
static void carry_prior(Prices& P, int nch) {
  u8 lt = P.ltype[P.nch - 1], ct = P.ctype[P.nch - 1], dt = P.dtype[P.nch - 1];
  memset(P.ltype, lt, nch); memset(P.ctype, ct, nch); memset(P.dtype, dt, nch); P.nch = nch;
}

// ---------------- context clustering ----------------
template <int D> struct HistT { u32 c[D]; u32 tot; };
static double g_hdr_nz = 8.0, g_hdr_base = 32.0;
static inline double hdr_est(int nz) { return g_hdr_base + g_hdr_nz * nz; }
template <int D> static double hist_cost(const HistT<D>& h) {
  double s = xlogx(h.tot); int nz = 0;
  for (int i = 0; i < D; i++) if (h.c[i]) { s -= xlogx(h.c[i]); nz++; }
  return s + hdr_est(nz);
}
// Greedy agglomerative merge of cl[first..first+n). Merges while it reduces cost or while more than
// maxcl remain. into[i] gets the surviving index for every input.
template <int D> static void merge_range(std::vector<HistT<D>>& cl, std::vector<double>& cc, int first, int n, int maxcl, std::vector<int>& into) {
  std::vector<double> Dm((size_t)n * n, 0.0);
  std::vector<char> alive(n, 1);
  auto pd = [&](int a, int b) {
    const HistT<D>& A = cl[first + a]; const HistT<D>& B = cl[first + b];
    double s = xlogx(A.tot + B.tot); int nz = 0;
    for (int i = 0; i < D; i++) { u32 v = A.c[i] + B.c[i]; if (v) { s -= xlogx(v); nz++; } }
    return s + hdr_est(nz) - cc[first + a] - cc[first + b];
  };
  for (int a = 0; a < n; a++) { into[first + a] = first + a; for (int b = a + 1; b < n; b++) Dm[(size_t)a * n + b] = pd(a, b); }
  int nalive = n;
  while (nalive > 1) {
    double best = 1e300; int ba = -1, bb = -1;
    for (int a = 0; a < n; a++) if (alive[a]) for (int b = a + 1; b < n; b++) if (alive[b]) {
      double d = Dm[(size_t)a * n + b]; if (d < best) { best = d; ba = a; bb = b; }
    }
    if (best >= 0 && nalive <= maxcl) break;
    HistT<D>& A = cl[first + ba]; const HistT<D>& B = cl[first + bb];
    for (int i = 0; i < D; i++) A.c[i] += B.c[i];
    A.tot += B.tot; cc[first + ba] = hist_cost(A); alive[bb] = 0; nalive--;
    for (int j = 0; j < n; j++) if (into[first + j] == first + bb) into[first + j] = first + ba;
    for (int j = 0; j < n; j++) if (alive[j] && j != ba) {
      int a = std::min(j, ba), b = std::max(j, ba); Dm[(size_t)a * n + b] = pd(a, b);
    }
  }
}
// returns number of clusters; fills map[ctx] and out[]; cost = estimated total bits
template <int D> static int cluster(const std::vector<HistT<D>>& h, int nctx, int maxcl, u8* map, std::vector<HistT<D>>& out, double& cost) {
  std::vector<int> owner(nctx, -1);
  std::vector<HistT<D>> cl; std::vector<double> cc;
  for (int i = 0; i < nctx; i++) if (h[i].tot) { owner[i] = (int)cl.size(); cl.push_back(h[i]); cc.push_back(hist_cost(h[i])); }
  int n = (int)cl.size();
  out.clear();
  if (n == 0) { memset(map, 0, nctx); HistT<D> z; memset(&z, 0, sizeof z); out.push_back(z); cost = 0; return 1; }
  std::vector<int> into(n);
  while (n > 96) {  // reduce in chunks, then compact survivors
    int before = n;
    for (int f = 0; f < n; f += 64) merge_range(cl, cc, f, std::min(64, n - f), 12, into);
    std::vector<int> newidx(n, -1); std::vector<HistT<D>> c2; std::vector<double> cc2;
    for (int i = 0; i < n; i++) if (into[i] == i) { newidx[i] = (int)c2.size(); c2.push_back(cl[i]); cc2.push_back(cc[i]); }
    for (int i = 0; i < nctx; i++) if (owner[i] >= 0) owner[i] = newidx[into[owner[i]]];
    cl.swap(c2); cc.swap(cc2); n = (int)cl.size(); into.assign(n, 0);
    if (n == before) break;
  }
  merge_range(cl, cc, 0, n, maxcl, into);
  std::vector<int> idx(n, -1); cost = 0;
  for (int a = 0; a < n; a++) if (into[a] == a) { idx[a] = (int)out.size(); out.push_back(cl[a]); cost += cc[a]; }
  int prev = 0, runs = 0;
  for (int i = 0; i < nctx; i++) {
    int v = owner[i] < 0 ? prev : idx[into[owner[i]]];
    if (i == 0 || v != prev) runs++;
    map[i] = (u8)v; prev = v;
  }
  int mb = 0; while ((1 << mb) < (int)out.size()) mb++;
  cost += (double)(mb + 6) * runs;
  return (int)out.size();
}

// ---------------- building models from a parse ----------------
struct LitInfo { std::vector<u8> v, p1, p2, mb, pos; std::vector<u16> ch; std::vector<u32> seqch; std::vector<u16> seqpos; };

static inline void rep_update(u32* r, u32 kind, u32 off) {
  if (kind == 0) return;
  if (kind == 1) { r[1] = r[0]; r[0] = off; return; }
  if (kind == 2) { r[2] = r[1]; r[1] = r[0]; r[0] = off; return; }
  r[3] = r[2]; r[2] = r[1]; r[1] = r[0]; r[0] = off;
}

// Walk a parse. With L: record every literal with its context bytes and chunk, and each sequence's
// (run-start chunk << 16 | match-start chunk).
static void replay(const u8* base, u32 bs, int lg, u32 PER, u32 start, u32 end, const u32* reps_in, const std::vector<Seq>& seqs, LitInfo* L, u32* reps_out) {
  u32 pos = start; u32 r[NUM_REPS]; memcpy(r, reps_in, sizeof r);
  if (L) { L->v.clear(); L->p1.clear(); L->p2.clear(); L->mb.clear(); L->pos.clear(); L->ch.clear(); L->seqch.clear(); L->seqpos.clear(); }
  auto lits = [&](u32 cnt) {
    if (L) for (u32 k = 0; k < cnt; k++, pos++) {
      const u8* q = base + pos;
      L->v.push_back(q[0]); L->p1.push_back(*(q - 1)); L->p2.push_back(*(q - 2)); L->mb.push_back(*(q - r[0]));
      L->pos.push_back((u8)(pos % PER)); L->ch.push_back((u16)((pos - bs) >> lg));
    } else pos += cnt;
  };
  for (const Seq& s : seqs) {
    u32 c0 = (pos - bs) >> lg, a0 = pos % PER;
    lits(s.ll);
    if (L) { L->seqch.push_back((c0 << 16) | ((pos - bs) >> lg)); L->seqpos.push_back((u16)(a0 | ((pos % PER) << 8))); }
    rep_update(r, s.kind, s.off); pos += s.ml;
  }
  lits(end - pos);
  if (reps_out) memcpy(reps_out, r, sizeof r);
}

static inline u32 off_symbol(const Seq& s, u32& nb, u32& ex) {  // for kind >= 3
  if (s.kind == 3) { u32 sym; off_code(s.off - 1, sym, nb, ex); return sym; }
  nb = 0; ex = 0; return OFF_SLOTS + (s.kind - 4);
}

// h holds ntypes * PER * sub histograms indexed [type][column][sub]. Clusters them as they are and
// with the position dimension summed out, keeps the cheaper, and fills map (full layout), freq tables
// and per-table symbol counts. Returns the number of tables.
template <int D> static int cluster_pos(const std::vector<HistT<D>>& h, int ntypes, int PER, int sub, int maxcl, u8* map,
                                        u16 (*tab)[D], u32* cnt, int bits, bool allow_pos) {
  const int nfull = ntypes * PER * sub, ncoll = ntypes * sub;
  std::vector<HistT<D>> hc(ncoll); memset(hc.data(), 0, sizeof(HistT<D>) * ncoll);
  for (int t = 0; t < ntypes; t++) for (int a = 0; a < PER; a++) for (int k = 0; k < sub; k++) {
    const HistT<D>& src = h[(t * PER + a) * sub + k]; HistT<D>& dst = hc[t * sub + k];
    if (!src.tot) continue;
    for (int i = 0; i < D; i++) dst.c[i] += src.c[i];
    dst.tot += src.tot;
  }
  std::vector<HistT<D>> of, oc; std::vector<u8> mf(nfull), mc(ncoll); double cf = 1e300, cc = 0;
  cluster(hc, ncoll, maxcl, mc.data(), oc, cc);
  if (allow_pos) cluster(h, nfull, maxcl, mf.data(), of, cf);
  const std::vector<HistT<D>>* o;
  if (cf < cc) { o = &of; memcpy(map, mf.data(), nfull); }
  else { o = &oc; for (int t = 0; t < ntypes; t++) for (int a = 0; a < PER; a++) for (int k = 0; k < sub; k++) map[(t * PER + a) * sub + k] = mc[t * sub + k]; }
  for (size_t c = 0; c < o->size(); c++) { normalize_freqs((*o)[c].c, D, tab[c], bits); if (g_coarse) coarsen_freqs((*o)[c].c, D, tab[c], bits); cnt[c] = (*o)[c].tot; }
  return (int)o->size();
}

template <int D> static std::vector<HistT<D>> zeros(int n) { std::vector<HistT<D>> h(n); memset(h.data(), 0, sizeof(HistT<D>) * n); return h; }

// Builds models for the block [bs, bs+bn) from the sequences covering [start, end) within it.
static void build_models(const u8* base, u32 bs, u32 bn, u32 start, u32 end, const u32* reps_in, const std::vector<Seq>& seqs,
                         const Params& prm, Models& M, LitInfo& L) {
  const int lg = prm.lgchunk; const int nch = (int)((bn + (1u << lg) - 1) >> lg);
  M.lgchunk = lg; M.nch = nch;
  const int PER = M.period;
  replay(base, bs, lg, (u32)PER, start, end, reps_in, seqs, &L, nullptr);
  const size_t nl = L.v.size(), ns = seqs.size();
  double cost;
  // --- chunk types
  {
    auto hl = zeros<256>(nch); auto hc = zeros<LL_SYMS + MLK_SYMS>(nch); auto hd = zeros<OFF_SYMS + ALIGN_SYMS>(nch);
    for (size_t i = 0; i < nl; i++) { auto& h = hl[L.ch[i]]; h.c[L.v[i]]++; h.tot++; }
    for (size_t i = 0; i < ns; i++) {
      const Seq& s = seqs[i]; u32 c0 = L.seqch[i] >> 16, c1 = L.seqch[i] & 0xffff; u32 sym, nb, ex;
      ll_code(s.ll, sym, nb, ex); hc[c0].c[sym]++; hc[c0].tot++;
      ml_code(s.ml - MIN_MATCH, sym, nb, ex); hc[c1].c[LL_SYMS + std::min<u32>(s.kind, 3) * MLC_SYMS + sym]++; hc[c1].tot++;
      if (s.kind >= 3) {
        sym = off_symbol(s, nb, ex); hd[c1].c[sym]++; hd[c1].tot++;
        if (nb >= 4) { hd[c1].c[OFF_SYMS + (ex & 15)]++; hd[c1].tot++; }
      }
    }
    std::vector<HistT<256>> o1; std::vector<HistT<LL_SYMS + MLK_SYMS>> o2; std::vector<HistT<OFF_SYMS + ALIGN_SYMS>> o3;
    M.nlt = cluster(hl, nch, prm.maxtypes, M.ltype, o1, cost);
    M.nct = cluster(hc, nch, prm.maxtypes, M.ctype, o2, cost);
    M.ndt = cluster(hd, nch, prm.maxtypes, M.dtype, o3, cost);
  }
  // --- literal mode per type: rank modes by a coarse clustering (the busiest contexts kept apart,
  //     the rest pooled), then cluster all (type, context) histograms jointly into shared tables.
  {
    std::vector<std::vector<u32>> byt(M.nlt);
    for (size_t i = 0; i < nl; i++) byt[M.ltype[L.ch[i]]].push_back((u32)i);
    const int RS = std::max(2, std::min(prm.rankslots, 64));
    std::vector<u32> tot; std::vector<u16> slot; std::vector<int> order; std::vector<HistT<256>> hs, outh; u8 smap[64];
    u32 nmap = 0;
    for (int t = 0; t < M.nlt; t++) {
      const std::vector<u32>& ix = byt[t];
      double bestc = 1e300; int bestmode = 0;
      for (int mode = 0; mode < NUM_LITMODES; mode++) {
        if (!((prm.modemask >> mode) & 1)) continue;
        const LitMode& lm = lit_mode(mode); const int nctx = lm.nctx;
        tot.assign(nctx, 0);
        for (u32 i : ix) tot[lm.A[L.p1[i]] + lm.B[L.mb[i]] + lm.C[L.p2[i]] + lm.D[L.pos[i]]]++;
        order.clear(); for (int c = 0; c < nctx; c++) if (tot[c]) order.push_back(c);
        int keep = std::min<int>((int)order.size(), RS - 1);
        if ((int)order.size() > keep) std::nth_element(order.begin(), order.begin() + keep, order.end(), [&](int a, int b) { return tot[a] > tot[b]; });
        slot.assign(nctx, (u16)keep);
        for (int k = 0; k < keep; k++) slot[order[k]] = (u16)k;
        hs.resize(keep + 1); memset(hs.data(), 0, sizeof(HistT<256>) * (keep + 1));
        for (u32 i : ix) {
          auto& hh = hs[slot[lm.A[L.p1[i]] + lm.B[L.mb[i]] + lm.C[L.p2[i]] + lm.D[L.pos[i]]]];
          hh.c[lm.sub ? (u8)(L.v[i] - L.mb[i]) : L.v[i]]++; hh.tot++;
        }
        cluster(hs, keep + 1, RS, smap, outh, cost);
        if (cost < bestc) { bestc = cost; bestmode = mode; }
      }
      M.lmode[t] = (u8)bestmode; M.lbase[t] = nmap; nmap += lit_mode(bestmode).nctx;
    }
    std::vector<HistT<256>> h(nmap); memset(h.data(), 0, sizeof(HistT<256>) * nmap);
    for (size_t i = 0; i < nl; i++) {
      u32 t = M.ltype[L.ch[i]]; const LitMode& lm = lit_mode(M.lmode[t]);
      auto& hh = h[M.lbase[t] + lm.A[L.p1[i]] + lm.B[L.mb[i]] + lm.C[L.p2[i]] + lm.D[L.pos[i]]];
      hh.c[lm.sub ? (u8)(L.v[i] - L.mb[i]) : L.v[i]]++; hh.tot++;
    }
    cluster(h, (int)nmap, prm.maxcl, M.litmap, outh, cost);
    M.ncl = (int)outh.size();
    for (int c = 0; c < M.ncl; c++) { normalize_freqs(outh[c].c, 256, M.lit[c], LIT_BITS); if (g_coarse) coarsen_freqs(outh[c].c, 256, M.lit[c], LIT_BITS); M.litn[c] = outh[c].tot; }
  }
  // --- command and distance tables
  auto hll = zeros<LL_SYMS>(M.nct * PER); auto hmlk = zeros<MLK_SYMS>(M.nct * PER * MLK_LLCTX);
  auto hoff = zeros<OFF_SYMS>(M.ndt * PER * OFF_LENCTX); auto hal = zeros<ALIGN_SYMS>(M.ndt);
  for (size_t i = 0; i < ns; i++) {
    const Seq& s = seqs[i]; u32 ct0 = M.ctype[L.seqch[i] >> 16], c1 = L.seqch[i] & 0xffff, ct1 = M.ctype[c1], dt1 = M.dtype[c1];
    u32 sym, nb, ex;
    const u32 a0 = L.seqpos[i] & 0xff, a1 = L.seqpos[i] >> 8;
    ll_code(s.ll, sym, nb, ex); { auto& hh = hll[ct0 * PER + a0]; hh.c[sym]++; hh.tot++; }
    ml_code(s.ml - MIN_MATCH, sym, nb, ex);
    { auto& hh = hmlk[(ct1 * PER + a1) * MLK_LLCTX + mlk_ctx(s.ll)]; hh.c[std::min<u32>(s.kind, 3) * MLC_SYMS + sym]++; hh.tot++; }
    if (s.kind >= 3) {
      sym = off_symbol(s, nb, ex);
      { auto& hh = hoff[(dt1 * PER + a1) * OFF_LENCTX + std::min<u32>(s.ml - MIN_MATCH, 3)]; hh.c[sym]++; hh.tot++; }
      if (nb >= 4) { hal[dt1].c[ex & 15]++; hal[dt1].tot++; }
    }
  }
  // Each of these is clustered twice, with and without the position-alignment part of the context,
  // and the cheaper result kept: alignment helps record-structured data and only dilutes text.
  M.nll = cluster_pos<LL_SYMS>(hll, M.nct, PER, 1, LL_MAXCL, M.llmap, M.ll, M.lln, LL_BITS, prm.posctx & 1);
  M.nmk = cluster_pos<MLK_SYMS>(hmlk, M.nct, PER, MLK_LLCTX, MLK_MAXCL, M.mlkmap, M.mlk, M.mlkn, MLK_BITS, prm.posctx & 2);
  M.nof = cluster_pos<OFF_SYMS>(hoff, M.ndt, PER, OFF_LENCTX, OFF_MAXCL, M.offmap, M.off, M.offn, OFF_BITS, prm.posctx & 4);
  { std::vector<HistT<ALIGN_SYMS>> o; M.nal = cluster(hal, M.ndt, AL_MAXCL, M.almap, o, cost);
    u64 na = 0; for (int c = 0; c < M.nal; c++) na += o[c].tot;
    // If modelling the low offset bits saves under 1.5%, store them raw: one less decode per offset.
    if (cost >= 4.0 * na * prm.align_gain) M.nal = 0;
    for (int c = 0; c < M.nal; c++) { normalize_freqs(o[c].c, ALIGN_SYMS, M.align[c], AL_BITS); M.aln[c] = o[c].tot; } }
}

// ---------------- optimal parse ----------------
struct Node { float cost; u32 ml, off, ll, kind; u32 rep[NUM_REPS]; };

// Lengths worth trying for a match of length `len`: every length up to `full`, then only the last
// length of each length-code bucket (the price is constant inside a bucket) and `len` itself.
static inline u32 next_len(u32 l, u32 len, u32 full) {
  if (l < full || l + 1 >= len) return l + 1;
  u32 v = l + 1 - MIN_MATCH;                 // code value of the next length
  if (v < 32) return l + 1;
  u32 n = ilog2(v); u32 half = 1u << (n - 1);
  u32 end = ((v >> (n - 1)) + 1) * half - 1; // last value in v's bucket
  u32 r = end + MIN_MATCH;
  return r < len ? r : len;
}

static inline float pll(const float* t, u32 v) { u32 s, nb, ex; ll_code(v, s, nb, ex); return t[s] + nb; }

// Parses the sub-block [bs+sb, bs+sb+bn) of the block starting at bs. carry_ll literals precede it
// unparsed; they are folded into the first sequence. Returns trailing literal count via carry_out.
static void parse(const u8* base, u32 bs, u32 sb, u32 bn, const u32* reps_in, u32 carry_ll, const MatchStore& ms,
                  const Prices& P, const Params& prm, std::vector<Node>& nodes, std::vector<u32>& mst, std::vector<Seq>& seqs,
                  u32* carry_out, u32* reps_out) {
  const u32 SEG = 1u << prm.lgseg;
  const u32 nice = prm.nice; const int lg = P.lgchunk;
  const u32 t2 = prm.far2, t3 = prm.far3; const u32 full = (u32)prm.fulllen;
  const u32 PER = (u32)P.period; std::vector<u8> colv;
  seqs.clear();
  u32 reps[NUM_REPS]; memcpy(reps, reps_in, sizeof reps);
  size_t mbase = 0;
  std::vector<Seq> rev;
  for (u32 ss = 0; ss < bn; ss += SEG) {
    const u32 sn = std::min(SEG, bn - ss);
    mst.resize(sn + 1);
    { size_t idx = mbase; for (u32 i = 0; i < sn; i++) { mst[i] = (u32)idx; idx += ms.cnt[ss + i]; } mst[sn] = (u32)idx; mbase = idx; }
    nodes.resize(sn + 1);
    for (u32 i = 0; i <= sn; i++) nodes[i].cost = INF;
    colv.resize(sn + 1); { u32 c = (bs + sb + ss) % PER; for (u32 i = 0; i <= sn; i++) { colv[i] = (u8)c; if (++c == PER) c = 0; } }
    const u8* col = colv.data();
    Node& n0 = nodes[0]; n0.cost = 0; n0.ml = 0; n0.off = 0; n0.ll = carry_ll; n0.kind = 0;
    memcpy(n0.rep, reps, sizeof n0.rep);
    auto relax = [&](u32 t, float c, const Node& N, u32 ml, u32 off, u32 kind) {
      Node& T = nodes[t];
      if (c < T.cost) {
        T.cost = c; T.ml = ml; T.off = off; T.kind = kind; T.ll = 0;
        memcpy(T.rep, N.rep, sizeof T.rep); rep_update(T.rep, kind, off);
      }
    };
    for (u32 i = 0; i < sn; i++) {
      const Node N = nodes[i];
      if (N.cost >= INF) continue;
      const u32 ch = std::min<u32>((sb + ss + i) >> lg, P.nch - 1);
      const u32 a = bs + sb + ss + i;
      const LitMode& lm = lit_mode(P.lmode[P.ltype[ch]]); const u8* lmap = P.litmap + P.lbase[P.ltype[ch]];
      const u8* llm = P.llmap + P.ctype[ch] * PER;
      const float* pl = P.ll[llm[(a - N.ll) % PER]];
      const float* pm = P.mlk[P.mlkmap[(P.ctype[ch] * PER + col[i]) * MLK_LLCTX + mlk_ctx(N.ll)]];
      const u8* omap = P.offmap + (P.dtype[ch] * PER + col[i]) * OFF_LENCTX;
      const float* pal = P.align[P.almap[P.dtype[ch]]];
      const u8* cur = base + a; u32 maxl = std::min<u32>(sn - i, MAX_MATCH);
      u32 rl[3] = {0, 0, 0}; u32 bestr = 0, bestk = 0;
      if (maxl >= 2) for (u32 k = 0; k < 3; k++) {
        if ((k >= 1 && N.rep[k] == N.rep[0]) || (k == 2 && N.rep[2] == N.rep[1])) continue;
        rl[k] = match_len(cur, cur - N.rep[k], maxl);
        if (rl[k] > bestr) { bestr = rl[k]; bestk = k; }
      }
      const u32 m0 = mst[i], nm = mst[i + 1] - m0;
      const u32* mo = ms.off.data() + m0; const u16* mln = ms.len.data() + m0;
      u32 bestm = nm ? std::min<u32>(mln[nm - 1], maxl) : 0;
      // cost of the next sequence's empty literal run, by the alignment of where this match ends
      const u8* colp = col + i;
#define LL0(l) (P.ll[llm[colp[l]]][0])
      const float base_c = N.cost;
      auto match_price = [&](u32 ml, u32 kind) { u32 s, nb, ex; ml_code(ml - MIN_MATCH, s, nb, ex); return pm[kind * MLC_SYMS + s] + nb; };
      auto off_price = [&](u32 off, u32 c) {
        u32 s, nb, ex; off_code(off - 1, s, nb, ex);
        return P.off[omap[c]][s] + (nb >= 4 ? (float)(nb - 4) + pal[ex & 15] : (float)nb);
      };
      if (bestr >= nice || bestm >= nice) {
        if (bestr >= bestm) relax(i + bestr, base_c + LL0(bestr) + match_price(bestr, bestk), N, bestr, N.rep[bestk], bestk);
        else relax(i + bestm, base_c + LL0(bestm) + match_price(bestm, 3) + off_price(mo[nm - 1], 3), N, bestm, mo[nm - 1], 3);
        i += std::max(bestr, bestm) - 1;
        continue;
      }
      {  // literal
        u8 mb = *(cur - N.rep[0]);
        u32 ctx = lm.A[cur[-1]] + lm.B[mb] + lm.C[cur[-2]] + lm.D[col[i]];
        u8 s = lm.sub ? (u8)(cur[0] - mb) : cur[0];
        float c = N.cost + P.lit[lmap[ctx]][s] + pll(pl, N.ll + 1) - pll(pl, N.ll);
        Node& T = nodes[i + 1];
        if (c < T.cost) { T = N; T.cost = c; T.ml = 0; T.ll = N.ll + 1; }
      }
      if (maxl < 2) continue;
      for (u32 k = 0; k < 3; k++) {
        for (u32 l = MIN_MATCH; l <= rl[k]; l = next_len(l, rl[k], full)) relax(i + l, base_c + LL0(l) + match_price(l, k), N, l, N.rep[k], k);
      }
      // special offsets: rep3 and rep0/rep1 +- 1..3
      const u16 c16 = (u16)load16(cur);
      for (u32 sp = 0; sp < NUM_SPECIAL; sp++) {
        u32 off;
        if (sp == 0) off = N.rep[3];
        else { int o = (int)N.rep[(sp - 1) / 6] + NEAR_DELTA[(sp - 1) % 6]; if (o < 1) continue; off = (u32)o; }
        if (off > a || off == N.rep[0] || off == N.rep[1] || off == N.rep[2]) continue;
        if (load16(cur - off) != c16) continue;
        u32 len = match_len(cur, cur - off, maxl);
        for (u32 l = MIN_MATCH; l <= len; l = next_len(l, len, full))
          relax(i + l, base_c + LL0(l) + match_price(l, 3) + P.off[omap[std::min<u32>(l - MIN_MATCH, 3)]][OFF_SLOTS + sp], N, l, off, 4 + sp);
      }
      u32 prev = MIN_MATCH - 1;
      for (u32 j = 0; j < nm; j++) {
        u32 lc = std::min<u32>(mln[j], maxl); u32 off = mo[j];
        if (lc <= prev) continue;
        if (off == N.rep[0] || off == N.rep[1] || off == N.rep[2]) { prev = lc; continue; }
        u32 lmin = off > t3 ? 4 : off > t2 ? 3 : 2;
        float op[4]; for (u32 c = 0; c < 4; c++) op[c] = off_price(off, c);
        for (u32 l = std::max(prev + 1, lmin); l <= lc; l = next_len(l, lc, full))
          relax(i + l, base_c + LL0(l) + match_price(l, 3) + op[std::min<u32>(l - MIN_MATCH, 3)], N, l, off, 3);
        prev = lc;
      }
    }
#undef LL0
    // backtrack this segment
    rev.clear(); u32 run = 0; u32 i = sn; u32 tail = 0; bool any = false;
    while (i > 0) {
      const Node& n = nodes[i];
      if (n.ml == 0) { run++; i--; }
      else {
        if (any) rev.back().ll = run; else { tail = run; any = true; }
        rev.push_back(Seq{0, n.ml, n.off, n.kind}); run = 0; i -= n.ml;
      }
    }
    if (any) {
      rev.back().ll = run + carry_ll; carry_ll = tail;
      seqs.insert(seqs.end(), rev.rbegin(), rev.rend());
    } else carry_ll += sn;
    memcpy(reps, nodes[sn].rep, sizeof reps);
  }
  *carry_out = carry_ll; memcpy(reps_out, reps, sizeof reps);
}

// ---------------- block emission ----------------
// One coded symbol: start (12 bits) | freq (12) | rANS state index (3) | precision (4)
typedef u32 Ev;
static double g_bits[8];  // lit, ll, mlk, off, align, len-extra, off-extra
static u64 g_kinds[20];
static u64 g_tables[8];
static u64 g_period[PMAX + 1];
static u64 g_hdr[8];  // header bits: type maps, lit map, lit tables, ll, mlk, off, align

static void cumul(const u16* f, int n, u16* cum) { u32 c = 0; for (int i = 0; i < n; i++) { cum[i] = (u16)c; c += f[i]; } }

static void emit_block(const u8* base, u32 bs, u32 bn, const u32* reps_in, const std::vector<Seq>& seqs,
                       const Models& M, std::vector<u8>& out, size_t* hdr_bytes) {
  static u16 clit[MAXCL][256], cll[LL_MAXCL][LL_SYMS], cmlk[MLK_MAXCL][MLK_SYMS], coff[OFF_MAXCL][OFF_SYMS], cal[AL_MAXCL][ALIGN_SYMS];
  for (int c = 0; c < M.ncl; c++) cumul(M.lit[c], 256, clit[c]);
  for (int c = 0; c < M.nll; c++) cumul(M.ll[c], LL_SYMS, cll[c]);
  for (int c = 0; c < M.nmk; c++) cumul(M.mlk[c], MLK_SYMS, cmlk[c]);
  for (int c = 0; c < M.nof; c++) cumul(M.off[c], OFF_SYMS, coff[c]);
  for (int c = 0; c < M.nal; c++) cumul(M.align[c], ALIGN_SYMS, cal[c]);
  std::vector<Ev> ev; ev.reserve(bn + 4 * seqs.size() + 16);
  BitWriter bits;
  auto add = [&](u16 start, u16 freq, u8 st, u8 bits) {
    if (!freq) { fprintf(stderr, "internal error: zero-frequency symbol\n"); abort(); }
    ev.push_back((u32)start | ((u32)freq << 12) | ((u32)st << 24) | ((u32)bits << 27));
  };
  u32 pos = bs, li = 0; u32 r[NUM_REPS]; memcpy(r, reps_in, sizeof r);
  const int lg = M.lgchunk; const LitMode* lm = &lit_mode(0); const u32 PER = (u32)M.period;
  const u32 end = bs + bn; u32 chunk_end = bs;
  const u8* lmap = nullptr; const u8* mmap = nullptr; const u8* omap = nullptr; const u8* llm = nullptr; u32 alc = 0;
  auto upd = [&]() {
    if (pos >= end) return;
    u32 ch = (pos - bs) >> lg; chunk_end = bs + ((ch + 1) << lg);
    lm = &lit_mode(M.lmode[M.ltype[ch]]); lmap = M.litmap + M.lbase[M.ltype[ch]]; llm = M.llmap + M.ctype[ch] * PER;
    mmap = M.mlkmap + M.ctype[ch] * PER * MLK_LLCTX; omap = M.offmap + M.dtype[ch] * PER * OFF_LENCTX; alc = M.almap[M.dtype[ch]];
  };
  auto lits = [&](u32 cnt) {
    while (cnt) {
      u32 n = std::min(cnt, chunk_end - pos);
      for (u32 k = 0; k < n; k++, pos++) {
        const u8* q = base + pos; u8 mb = *(q - r[0]); u32 cl = lmap[lm->A[*(q - 1)] + lm->B[mb] + lm->C[*(q - 2)] + lm->D[pos % PER]];
        u8 s = lm->sub ? (u8)(base[pos] - mb) : base[pos];
        add(clit[cl][s], M.lit[cl][s], (u8)(li & 1), LIT_BITS); li++; g_bits[0] += FPRICE[(u32)M.lit[cl][s] << (MAX_BITS - LIT_BITS)];
      }
      cnt -= n; if (pos >= chunk_end) upd();
    }
  };
  for (const Seq& s : seqs) {
    u32 sym, nb, ex;
    if (pos >= chunk_end) upd();
    const u32 llc = llm[pos % PER];
    ll_code(s.ll, sym, nb, ex); add(cll[llc][sym], M.ll[llc][sym], 2, LL_BITS); bits.put(ex, nb); g_bits[1] += FPRICE[(u32)M.ll[llc][sym] << (MAX_BITS - LL_BITS)]; g_bits[5] += nb;
    lits(s.ll);
    ml_code(s.ml - MIN_MATCH, sym, nb, ex);
    const u32 mcol = pos % PER;
    u32 c = mmap[mcol * MLK_LLCTX + mlk_ctx(s.ll)]; u32 ms = std::min<u32>(s.kind, 3) * MLC_SYMS + sym;
    add(cmlk[c][ms], M.mlk[c][ms], 3, MLK_BITS); bits.put(ex, nb); g_bits[2] += FPRICE[(u32)M.mlk[c][ms] << (MAX_BITS - MLK_BITS)]; g_bits[5] += nb;
    g_kinds[s.kind]++;
    if (s.kind >= 3) {
      sym = off_symbol(s, nb, ex);
      u32 oc = omap[mcol * OFF_LENCTX + std::min<u32>(s.ml - MIN_MATCH, 3)];
      add(coff[oc][sym], M.off[oc][sym], 4, OFF_BITS); g_bits[3] += FPRICE[(u32)M.off[oc][sym] << (MAX_BITS - OFF_BITS)]; g_bits[6] += nb >= 4 ? nb - 4 : nb;
      if (nb >= 4 && M.nal) { bits.put(ex >> 4, nb - 4); add(cal[alc][ex & 15], M.align[alc][ex & 15], 3, AL_BITS); g_bits[4] += FPRICE[(u32)M.align[alc][ex & 15] << (MAX_BITS - AL_BITS)]; }
      else if (nb >= 4) { bits.put(ex, nb); g_bits[4] += 4; }
      else bits.put(ex, nb);
    }
    rep_update(r, s.kind, s.off);
    pos += s.ml;
  }
  if (pos >= chunk_end) upd();
  lits(end - pos);
  bits.flush();
  // rANS encode backwards
  std::vector<u16> w(ev.size() + 2 * NUM_STREAMS + 4);
  u16* ptr = w.data() + w.size();
  u32 x[NUM_STREAMS]; for (int i = 0; i < NUM_STREAMS; i++) x[i] = RANS_L;
  for (size_t i = ev.size(); i-- > 0;) {
    const Ev e = ev[i]; const u32 st = (e >> 24) & 7, eb = e >> 27, fr = (e >> 12) & 0xfff; u32 xs = x[st];
    u64 xmax = (u64)((RANS_L >> eb) << 16) * fr;
    if (xs >= xmax) { *--ptr = (u16)xs; xs >>= 16; }
    x[st] = ((xs / fr) << eb) + (xs % fr) + (e & 0xfff);
  }
  for (int i = NUM_STREAMS - 1; i >= 0; i--) { *--ptr = (u16)(x[i] >> 16); *--ptr = (u16)x[i]; }
  size_t rb = (w.data() + w.size() - ptr) * 2;
  // header
  RcEnc hw; MapModel mm[8]; TabModel tml, tmll, tmmk, tmof, tmal;
  auto hbits = [&]() { return (u64)hw.out.size() * 8; }; u64 h0 = 0;
#define HSTAT(k) do { u64 h1 = hbits(); g_hdr[k] += h1 - h0; h0 = h1; } while (0)
  hw.direct(lg, 5); hw.direct(PER, 6);
  hw.direct(M.nlt - 1, 4); write_map(hw, mm[0], M.ltype, M.nch, M.nlt);
  hw.direct(M.nct - 1, 4); write_map(hw, mm[1], M.ctype, M.nch, M.nct);
  hw.direct(M.ndt - 1, 4); write_map(hw, mm[2], M.dtype, M.nch, M.ndt);
  HSTAT(0);
  u32 nmap = 0; for (int t = 0; t < M.nlt; t++) { hw.direct(M.lmode[t], 4); nmap += lit_mode(M.lmode[t]).nctx; }
  hw.direct(M.ncl - 1, 6); write_map(hw, mm[3], M.litmap, (int)nmap, M.ncl);
  HSTAT(1);
  for (int c = 0; c < M.ncl; c++) write_table(hw, tml, M.lit[c], 256, LIT_BITS);
  HSTAT(2);
  hw.direct(M.nll - 1, 3); write_map(hw, mm[4], M.llmap, M.nct * PER, M.nll);
  for (int c = 0; c < M.nll; c++) write_table(hw, tmll, M.ll[c], LL_SYMS, LL_BITS);
  HSTAT(3);
  hw.direct(M.nmk - 1, 4); write_map(hw, mm[5], M.mlkmap, M.nct * PER * MLK_LLCTX, M.nmk);
  for (int c = 0; c < M.nmk; c++) write_table(hw, tmmk, M.mlk[c], MLK_SYMS, MLK_BITS);
  HSTAT(4);
  hw.direct(M.nof - 1, 4); write_map(hw, mm[6], M.offmap, M.ndt * PER * OFF_LENCTX, M.nof);
  for (int c = 0; c < M.nof; c++) write_table(hw, tmof, M.off[c], OFF_SYMS, OFF_BITS);
  HSTAT(5);
  hw.direct(M.nal, 3); if (M.nal) write_map(hw, mm[7], M.almap, M.ndt, M.nal);
  for (int c = 0; c < M.nal; c++) write_table(hw, tmal, M.align[c], ALIGN_SYMS, AL_BITS);
  hw.finish();
  HSTAT(6);
#undef HSTAT
  g_tables[0] += M.ncl; g_tables[1] += M.nll; g_tables[2] += M.nmk; g_tables[3] += M.nof; g_tables[4] += M.nlt; g_tables[5] += M.nct; g_tables[6] += M.ndt; g_tables[7]++;
  put_varint(out, bn); out.push_back(1);
  put_varint(out, seqs.size()); put_varint(out, hw.out.size()); put_varint(out, rb); put_varint(out, bits.buf.size());
  out.insert(out.end(), hw.out.begin(), hw.out.end());
  out.insert(out.end(), (const u8*)ptr, (const u8*)ptr + rb);
  out.insert(out.end(), bits.buf.begin(), bits.buf.end());
  if (hdr_bytes) *hdr_bytes += hw.out.size();
}

// Column period of a block: the P in 1..PMAX for which byte values depend most on (position mod P).
// 4 unless some other period is clearly better (fixed-width records, 3-byte pixels, ...).
static int detect_period(const u8* p, u32 n, u32 filepos) {
  static u32 h[PMAX * 256];
  double ent[PMAX + 1];
  for (int P = 1; P <= PMAX; P++) {
    memset(h, 0, sizeof(u32) * P * 256);
    u32 c = filepos % P;
    for (u32 i = 0; i < n; i++) { h[c * 256 + p[i]]++; if (++c == (u32)P) c = 0; }
    double e = 0;
    for (int k = 0; k < P; k++) {
      u32 tot = 0; const u32* hh = h + k * 256;
      for (int s = 0; s < 256; s++) tot += hh[s];
      e += xlogx(tot); for (int s = 0; s < 256; s++) e -= xlogx(hh[s]);
    }
    ent[P] = e + 8.0 * 64 * P;  // small charge per column against fitting noise
  }
  int best = 4;
  for (int P = 1; P <= PMAX; P++) if (ent[P] < ent[best]) best = P;
  if (best == 4 || ent[best] > 0.97 * ent[4]) return 4;
  for (int P = 1; P <= PMAX; P++) if (ent[P] <= ent[best] * 1.005) return P;
  return best;
}

// Streaming compressor. Holds one window of history plus the block being coded and some lookahead,
// not the whole file.
bool compress(FILE* fin, u64 n, FILE* fout, const Params& prm, u64* out_size) {
  if (n >= 0xfff00000ull) { fprintf(stderr, "inputs of 4 GB or more are not supported\n"); return false; }
  einit(); g_coarse = prm.coarse; g_smooth = prm.smooth; g_unseen_cap = prm.unseen_cap; g_hdr_nz = prm.hdr_nz; g_hdr_base = prm.hdr_base;
  MF mf(n, prm);
  const u32 B = 1u << prm.lgblock;
  const u64 LOOK = (u64)std::max(prm.nice, 1 << 16) + 64;
  const u64 HIST = (u64)mf.W + 64;                       // history that must stay addressable
  const u64 cap = std::min<u64>(n, HIST + B + mf.W / 4 + LOOK);
  std::vector<u8> buf(cap + 2 * PAD + 96, 0);
  u8* data = buf.data() + PAD;
  u64 win_start = 0, loaded = 0;                         // buf holds file positions [win_start, loaded)
  const u8* base = data;
  auto fill = [&]() {
    size_t want = (size_t)std::min<u64>(n - loaded, cap - (loaded - win_start));
    size_t got = want ? fread(data + (loaded - win_start), 1, want, fin) : 0;
    loaded += got; memset(data + (loaded - win_start), 0, 64);
    mf.b = base; mf.n = (u32)loaded;
    return got == want;
  };
  if (!fill()) return false;
  std::vector<u8> out = {'G', 'T', '0', 8};
  put_varint(out, n);
  out.push_back((u8)ilog2(mf.W)); out.push_back((u8)prm.lgblock);
  u64 total_out = 0;
  Prices* pr = new Prices; default_prices(*pr, prm.lgchunk);
  Models* M = new Models;
  u32 reps[NUM_REPS] = {1, 2, 3, 4};
  const u32 SB = std::min(B, 1u << prm.lgsub);
  MatchStore ms; std::vector<u32> mst; std::vector<Match> tmp(prm.nice + 8);
  std::vector<Node> nodes; std::vector<Seq> seqs, sseqs; LitInfo L;
  int modestat[NUM_LITMODES] = {0}; size_t hdr_bytes = 0, nseq_tot = 0, stored = 0;
  const u32 t2 = prm.far2, t3 = prm.far3;
  bool ok = true;
  for (u64 bs64 = 0; bs64 < n && ok; bs64 += B) {
    const u32 bs = (u32)bs64; const u32 bn = (u32)std::min<u64>(B, n - bs64);
    if (loaded < std::min<u64>(n, bs64 + bn + LOOK)) {  // need more input: slide the window, then read
      u64 ns = bs64 > HIST ? bs64 - HIST : 0;
      if (ns > win_start) {
        memmove(data, data + (ns - win_start), (size_t)(loaded - ns));
        win_start = ns; base = data - win_start;
      }
      if (!fill()) { ok = false; break; }
    }
    const int PER = prm.period ? std::min(std::max(prm.period, 1), (int)PMAX) : detect_period(base + bs, std::min(bn, 1u << 18), bs);
    if (PER != pr->period) { default_prices(*pr, prm.lgchunk); pr->period = PER; }
    M->period = PER; g_period[PER]++;
    carry_prior(*pr, (int)((bn + (1u << prm.lgchunk) - 1) >> prm.lgchunk));
    seqs.clear();
    u32 sreps[NUM_REPS]; memcpy(sreps, reps, sizeof sreps); u32 carry = 0;
    for (u32 sb = 0; sb < bn; sb += SB) {
      const u32 sn = std::min(SB, bn - sb);
      double t0 = now();
      ms.cnt.assign(sn, 0); ms.off.clear(); ms.len.clear(); ms.off.reserve((size_t)sn * 6); ms.len.reserve((size_t)sn * 6);
      for (u32 i = 0; i < sn; i++) {
        u32 cnt = mf.find(bs + sb + i, tmp.data(), true);
        u32 kept = 0;
        for (u32 j = 0; j < cnt && kept < 255; j++) {
          u32 l = tmp[j].len, o = tmp[j].off;
          if ((l == 2 && o > t2) || (l == 3 && o > t3)) continue;
          ms.off.push_back(o); ms.len.push_back((u16)l); kept++;
        }
        ms.cnt[i] = (u8)kept;
        if (cnt && tmp[cnt - 1].len >= (u32)prm.nice) {
          u32 Lm = std::min(tmp[cnt - 1].len, sn - i);
          for (u32 j = 1; j < Lm; j++) mf.find(bs + sb + i + j, tmp.data(), false);
          i += Lm - 1;
        }
      }
      if (prm.verbose > 1) fprintf(stderr, "sub-block %u+%u: %u bytes, %zu matches (%.2f/pos)\n", bs, sb, sn, ms.off.size(), (double)ms.off.size() / sn);
      g_time[0] += now() - t0;
      int passes = prm.passes;
      u32 creps[NUM_REPS], cout = 0;
      const u32 mstart = bs + sb - carry;
      for (int it = 0; it < passes; it++) {
        double t0 = now();
        parse(base, bs, sb, sn, sreps, carry, ms, *pr, prm, nodes, mst, sseqs, &cout, creps);
        g_time[1] += now() - t0; t0 = now();
        if (SB < bn || it < passes - 1) {
          build_models(base, bs, bn, mstart, bs + sb + sn, sreps, sseqs, prm, *M, L);
          prices_from_models(*M, *pr);
        }
        g_time[2] += now() - t0;
      }
      seqs.insert(seqs.end(), sseqs.begin(), sseqs.end());
      carry = cout; memcpy(sreps, creps, sizeof sreps);
    }
    double t1 = now();
    build_models(base, bs, bn, bs, bs + bn, reps, seqs, prm, *M, L);
    for (int t = 0; t < M->nlt; t++) modestat[M->lmode[t]]++;
    prices_from_models(*M, *pr);
    g_time[2] += now() - t1; t1 = now();
    size_t before = out.size();
    emit_block(base, bs, bn, reps, seqs, *M, out, &hdr_bytes);
    g_time[3] += now() - t1;
    if (out.size() - before >= bn + 6) {
      out.resize(before); put_varint(out, bn); out.push_back(0);
      out.insert(out.end(), base + bs, base + bs + bn); stored++;
    } else {
      memcpy(reps, sreps, sizeof reps);
      nseq_tot += seqs.size();
    }
    if (fwrite(out.data(), 1, out.size(), fout) != out.size()) ok = false;
    total_out += out.size(); out.clear();
  }
  if (n == 0) { if (fwrite(out.data(), 1, out.size(), fout) != out.size()) ok = false; total_out += out.size(); }
  if (out_size) *out_size = total_out;
  if (prm.verbose) {
    fprintf(stderr, "time: matchfind %.2f parse %.2f models %.2f emit %.2f\n", g_time[0], g_time[1], g_time[2], g_time[3]);
    fprintf(stderr, "periods:"); for (int i = 1; i <= PMAX; i++) if (g_period[i]) fprintf(stderr, " %d x%llu", i, (unsigned long long)g_period[i]); fprintf(stderr, "\n");
    fprintf(stderr, "litmodes:"); for (int i = 0; i < NUM_LITMODES; i++) fprintf(stderr, " %d", modestat[i]);
    fprintf(stderr, "  stored %zu  hdr %zu B  seqs %zu\n", stored, hdr_bytes, nseq_tot);
    double nb = (double)std::max<u64>(1, g_tables[7]);
    fprintf(stderr, "avg/block: lit tables %.1f ll %.1f mlk %.1f off %.1f | types lit %.1f cmd %.1f dist %.1f\n", g_tables[0] / nb, g_tables[1] / nb, g_tables[2] / nb, g_tables[3] / nb, g_tables[4] / nb, g_tables[5] / nb, g_tables[6] / nb);
    fprintf(stderr, "header bytes: type-maps %llu lit-map %llu lit-tables %llu ll %llu mlk %llu off %llu align %llu\n", (unsigned long long)g_hdr[0] / 8, (unsigned long long)g_hdr[1] / 8, (unsigned long long)g_hdr[2] / 8, (unsigned long long)g_hdr[3] / 8, (unsigned long long)g_hdr[4] / 8, (unsigned long long)g_hdr[5] / 8, (unsigned long long)g_hdr[6] / 8);
    fprintf(stderr, "bytes: lit %.0f ll %.0f mlk %.0f off %.0f align %.0f len-extra %.0f off-extra %.0f\n", g_bits[0] / 8, g_bits[1] / 8, g_bits[2] / 8, g_bits[3] / 8, g_bits[4] / 8, g_bits[5] / 8, g_bits[6] / 8);
    fprintf(stderr, "kinds: rep0 %llu rep1 %llu rep2 %llu new %llu rep3 %llu near0", (unsigned long long)g_kinds[0], (unsigned long long)g_kinds[1], (unsigned long long)g_kinds[2], (unsigned long long)g_kinds[3], (unsigned long long)g_kinds[4]);
    u64 a = 0, b = 0; for (int i = 5; i < 11; i++) a += g_kinds[i]; for (int i = 11; i < 17; i++) b += g_kinds[i];
    fprintf(stderr, " %llu near1 %llu\n", (unsigned long long)a, (unsigned long long)b);
  }
  delete pr; delete M;
  return ok;
}

}  // namespace gt
