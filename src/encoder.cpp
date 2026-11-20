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
  MF(const u8* base, u32 n_, const Params& P) {
    b = base; n = n_; nice = P.nice; depth = P.depth; hb = P.hashbits;
    u32 w = 1u << P.lgwin; while (w / 2 >= n && w > 65536) w /= 2;
    W = w; wmask = W - 1;
    head.assign((size_t)1 << hb, NIL); son.resize(2 * (size_t)std::min<u64>(W, (u64)n + 1));
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
  int lgchunk, nch, nlt, nct, ndt;
  u8 ltype[MAXCH], ctype[MAXCH], dtype[MAXCH];
  int litmode, ncl, nll, nmk, nof, nal;
  u8 litmap[NT * MAXCTX], llmap[NT], mlkmap[NT * MLK_NCTX], offmap[NT * 4], almap[NT];
  u16 lit[MAXCL][256], ll[LL_MAXCL][LL_SYMS], mlk[MLK_MAXCL][MLK_SYMS], off[OFF_MAXCL][OFF_SYMS], align[AL_MAXCL][ALIGN_SYMS];
};
struct Prices {
  int lgchunk, nch, litmode;
  u8 ltype[MAXCH], ctype[MAXCH], dtype[MAXCH];
  u8 litmap[NT * MAXCTX], llmap[NT], mlkmap[NT * MLK_NCTX], offmap[NT * 4], almap[NT];
  float lit[MAXCL][256], ll[LL_MAXCL][LL_SYMS], mlk[MLK_MAXCL][MLK_SYMS], off[OFF_MAXCL][OFF_SYMS], align[AL_MAXCL][ALIGN_SYMS];
};
static float FPRICE[SCALE + 1];
static double XLOGX[65536];
static bool g_einit = false;
static void einit() {
  if (g_einit) return;
  FPRICE[0] = 13.0f; for (int i = 1; i <= SCALE; i++) FPRICE[i] = (float)log2((double)SCALE / i);
  XLOGX[0] = 0; for (int i = 1; i < 65536; i++) XLOGX[i] = i * log2((double)i);
  g_einit = true;
}
static inline double xlogx(u32 x) { return x < 65536 ? XLOGX[x] : x * log2((double)x); }

static void default_prices(Prices& P, int lgchunk) {
  memset(&P, 0, sizeof P);
  P.lgchunk = lgchunk; P.nch = 1;
  for (int s = 0; s < 256; s++) P.lit[0][s] = 6.5f;
  for (int s = 0; s < LL_SYMS; s++) P.ll[0][s] = 2.0f + 0.25f * s;
  for (int s = 0; s < MLK_SYMS; s++) {
    int k = s / MLC_SYMS, mc = s % MLC_SYMS;
    P.mlk[0][s] = (k == 0 ? 3.0f : k == 3 ? 4.0f : 6.0f) + 0.1f * std::min(mc, 32);
  }
  for (int s = 0; s < OFF_SYMS; s++) P.off[0][s] = s < OFF_SLOTS ? 5.5f : 8.0f;
  for (int s = 0; s < ALIGN_SYMS; s++) P.align[0][s] = 4.0f;
}
static void prices_from_models(const Models& M, Prices& P) {
  P.lgchunk = M.lgchunk; P.nch = M.nch; P.litmode = M.litmode;
  memcpy(P.ltype, M.ltype, sizeof P.ltype); memcpy(P.ctype, M.ctype, sizeof P.ctype); memcpy(P.dtype, M.dtype, sizeof P.dtype);
  memcpy(P.litmap, M.litmap, sizeof P.litmap); memcpy(P.llmap, M.llmap, sizeof P.llmap);
  memcpy(P.mlkmap, M.mlkmap, sizeof P.mlkmap); memcpy(P.offmap, M.offmap, sizeof P.offmap); memcpy(P.almap, M.almap, sizeof P.almap);
  for (int c = 0; c < M.ncl; c++) for (int s = 0; s < 256; s++) P.lit[c][s] = FPRICE[M.lit[c][s]];
  for (int c = 0; c < M.nll; c++) for (int s = 0; s < LL_SYMS; s++) P.ll[c][s] = FPRICE[M.ll[c][s]];
  for (int c = 0; c < M.nmk; c++) for (int s = 0; s < MLK_SYMS; s++) P.mlk[c][s] = FPRICE[M.mlk[c][s]];
  for (int c = 0; c < M.nof; c++) for (int s = 0; s < OFF_SYMS; s++) P.off[c][s] = FPRICE[M.off[c][s]];
  for (int c = 0; c < M.nal; c++) for (int s = 0; s < ALIGN_SYMS; s++) P.align[c][s] = FPRICE[M.align[c][s]];
}
// Reuse the previous block's statistics as the starting guess for a new block: every chunk of the
// new block takes the types of the previous block's last chunk.
static void carry_prior(Prices& P, int nch) {
  u8 lt = P.ltype[P.nch - 1], ct = P.ctype[P.nch - 1], dt = P.dtype[P.nch - 1];
  memset(P.ltype, lt, nch); memset(P.ctype, ct, nch); memset(P.dtype, dt, nch); P.nch = nch;
}

// ---------------- context clustering ----------------
template <int D> struct HistT { u32 c[D]; u32 tot; };
static inline double hdr_est(int nz) { return 32 + 8.0 * nz; }
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
struct LitInfo { std::vector<u8> v, p1, p2, mb, pos; std::vector<u16> ch; std::vector<u32> seqch; };

static inline void rep_update(u32* r, u32 kind, u32 off) {
  if (kind == 0) return;
  if (kind == 1) { r[1] = r[0]; r[0] = off; return; }
  if (kind == 2) { r[2] = r[1]; r[1] = r[0]; r[0] = off; return; }
  r[3] = r[2]; r[2] = r[1]; r[1] = r[0]; r[0] = off;
}

// Walk a parse. With L: record every literal with its context bytes and chunk, and each sequence's
// (run-start chunk << 16 | match-start chunk).
static void replay(const u8* base, u32 bs, int lg, u32 start, u32 end, const u32* reps_in, const std::vector<Seq>& seqs, LitInfo* L, u32* reps_out) {
  u32 pos = start; u32 r[NUM_REPS]; memcpy(r, reps_in, sizeof r);
  if (L) { L->v.clear(); L->p1.clear(); L->p2.clear(); L->mb.clear(); L->pos.clear(); L->ch.clear(); L->seqch.clear(); }
  auto lits = [&](u32 cnt) {
    if (L) for (u32 k = 0; k < cnt; k++, pos++) {
      const u8* q = base + pos;
      L->v.push_back(q[0]); L->p1.push_back(*(q - 1)); L->p2.push_back(*(q - 2)); L->mb.push_back(*(q - r[0]));
      L->pos.push_back((u8)(pos & 3)); L->ch.push_back((u16)((pos - bs) >> lg));
    } else pos += cnt;
  };
  for (const Seq& s : seqs) {
    u32 c0 = (pos - bs) >> lg;
    lits(s.ll);
    if (L) L->seqch.push_back((c0 << 16) | ((pos - bs) >> lg));
    rep_update(r, s.kind, s.off); pos += s.ml;
  }
  lits(end - pos);
  if (reps_out) memcpy(reps_out, r, sizeof r);
}

static inline u32 off_symbol(const Seq& s, u32& nb, u32& ex) {  // for kind >= 3
  if (s.kind == 3) { u32 sym; off_code(s.off - 1, sym, nb, ex); return sym; }
  nb = 0; ex = 0; return OFF_SLOTS + (s.kind - 4);
}

template <int D> static std::vector<HistT<D>> zeros(int n) { std::vector<HistT<D>> h(n); memset(h.data(), 0, sizeof(HistT<D>) * n); return h; }

// Builds models for the block [bs, bs+bn) from the sequences covering [start, end) within it.
static void build_models(const u8* base, u32 bs, u32 bn, u32 start, u32 end, const u32* reps_in, const std::vector<Seq>& seqs,
                         const Params& prm, Models& M, LitInfo& L, int* modestat, bool fullmodes) {
  const int lg = prm.lgchunk; const int nch = (int)((bn + (1u << lg) - 1) >> lg);
  M.lgchunk = lg; M.nch = nch;
  replay(base, bs, lg, start, end, reps_in, seqs, &L, nullptr);
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
  // --- literal mode: rank modes with a cheap adaptive-coding cost estimate, cluster only the best
  //     few to pick one (ignoring types), then do the typed context clustering for the winner.
  {
    std::vector<HistT<256>> h, outh, besth; static u8 map[NT * MAXCTX], bestmap[MAXCTX];
    auto fill = [&](int mode) {
      int nctx = LIT_NCTX[mode]; bool sub = lit_sub(mode);
      h.resize(nctx); memset(h.data(), 0, sizeof(HistT<256>) * nctx);
      for (size_t i = 0; i < nl; i++) {
        auto& hh = h[lit_ctx(mode, L.pos[i], L.p1[i], L.p2[i], L.mb[i])];
        hh.c[sub ? (u8)(L.v[i] - L.mb[i]) : L.v[i]]++; hh.tot++;
      }
    };
    int cand[NUM_LITMODES]; int ncand = 0;
    if (!fullmodes) cand[ncand++] = M.litmode;
    else if (prm.modetop >= NUM_LITMODES) { for (int m = 0; m < NUM_LITMODES; m++) if ((prm.modemask >> m) & 1) cand[ncand++] = m; }
    else {
      std::pair<double, int> pc[NUM_LITMODES]; int np = 0;
      const double al = prm.alpha;
      for (int mode = 0; mode < NUM_LITMODES; mode++) {
        if (!((prm.modemask >> mode) & 1)) continue;
        fill(mode);
        int nctx = LIT_NCTX[mode]; double c = 0;
        for (int x = 0; x < nctx; x++) if (h[x].tot) {
          const HistT<256>& hh = h[x];
          c += lgamma(hh.tot + 256 * al) - lgamma(256 * al);
          for (int k = 0; k < 256; k++) if (hh.c[k]) c -= lgamma(hh.c[k] + al) - lgamma(al);
        }
        pc[np++] = {c, mode};
      }
      std::sort(pc, pc + np);
      for (int k = 0; k < np && k < prm.modetop; k++) cand[ncand++] = pc[k].second;
    }
    double bestc = 1e300; int bestmode = cand[0];
    for (int k = 0; k < ncand; k++) {
      int mode = cand[k]; fill(mode);
      cluster(h, LIT_NCTX[mode], prm.maxcl, map, outh, cost);
      if (cost < bestc) { bestc = cost; bestmode = mode; besth = outh; memcpy(bestmap, map, LIT_NCTX[mode]); }
    }
    M.litmode = bestmode; if (modestat) modestat[bestmode]++;
    int nctx = LIT_NCTX[bestmode]; bool sub = lit_sub(bestmode);
    if (M.nlt == 1) { memcpy(M.litmap, bestmap, nctx); outh = besth; }
    else {
      h.resize((size_t)M.nlt * nctx); memset(h.data(), 0, sizeof(HistT<256>) * h.size());
      for (size_t i = 0; i < nl; i++) {
        auto& hh = h[(size_t)M.ltype[L.ch[i]] * nctx + lit_ctx(bestmode, L.pos[i], L.p1[i], L.p2[i], L.mb[i])];
        hh.c[sub ? (u8)(L.v[i] - L.mb[i]) : L.v[i]]++; hh.tot++;
      }
      cluster(h, M.nlt * nctx, prm.maxcl, M.litmap, outh, cost);
    }
    M.ncl = (int)outh.size();
    for (int c = 0; c < M.ncl; c++) normalize_freqs(outh[c].c, 256, M.lit[c]);
  }
  // --- command and distance tables
  auto hll = zeros<LL_SYMS>(M.nct); auto hmlk = zeros<MLK_SYMS>(M.nct * MLK_NCTX);
  auto hoff = zeros<OFF_SYMS>(M.ndt * 4); auto hal = zeros<ALIGN_SYMS>(M.ndt);
  for (size_t i = 0; i < ns; i++) {
    const Seq& s = seqs[i]; u32 ct0 = M.ctype[L.seqch[i] >> 16], c1 = L.seqch[i] & 0xffff, ct1 = M.ctype[c1], dt1 = M.dtype[c1];
    u32 sym, nb, ex;
    ll_code(s.ll, sym, nb, ex); hll[ct0].c[sym]++; hll[ct0].tot++;
    ml_code(s.ml - MIN_MATCH, sym, nb, ex);
    { auto& hh = hmlk[ct1 * MLK_NCTX + mlk_ctx(s.ll)]; hh.c[std::min<u32>(s.kind, 3) * MLC_SYMS + sym]++; hh.tot++; }
    if (s.kind >= 3) {
      sym = off_symbol(s, nb, ex);
      { auto& hh = hoff[dt1 * 4 + std::min<u32>(s.ml - MIN_MATCH, 3)]; hh.c[sym]++; hh.tot++; }
      if (nb >= 4) { hal[dt1].c[ex & 15]++; hal[dt1].tot++; }
    }
  }
  { std::vector<HistT<LL_SYMS>> o; M.nll = cluster(hll, M.nct, LL_MAXCL, M.llmap, o, cost);
    for (int c = 0; c < M.nll; c++) normalize_freqs(o[c].c, LL_SYMS, M.ll[c]); }
  { std::vector<HistT<MLK_SYMS>> o; M.nmk = cluster(hmlk, M.nct * MLK_NCTX, MLK_MAXCL, M.mlkmap, o, cost);
    for (int c = 0; c < M.nmk; c++) normalize_freqs(o[c].c, MLK_SYMS, M.mlk[c]); }
  { std::vector<HistT<OFF_SYMS>> o; M.nof = cluster(hoff, M.ndt * 4, OFF_MAXCL, M.offmap, o, cost);
    for (int c = 0; c < M.nof; c++) normalize_freqs(o[c].c, OFF_SYMS, M.off[c]); }
  { std::vector<HistT<ALIGN_SYMS>> o; M.nal = cluster(hal, M.ndt, AL_MAXCL, M.almap, o, cost);
    for (int c = 0; c < M.nal; c++) normalize_freqs(o[c].c, ALIGN_SYMS, M.align[c]); }
}

// ---------------- optimal parse ----------------
struct Node { float cost; u32 ml, off, ll, kind; u32 rep[NUM_REPS]; };

static inline float pll(const float* t, u32 v) { u32 s, nb, ex; ll_code(v, s, nb, ex); return t[s] + nb; }

// Parses the sub-block [bs+sb, bs+sb+bn) of the block starting at bs. carry_ll literals precede it
// unparsed; they are folded into the first sequence. Returns trailing literal count via carry_out.
static void parse(const u8* base, u32 bs, u32 sb, u32 bn, const u32* reps_in, u32 carry_ll, const MatchStore& ms,
                  const Prices& P, const Params& prm, std::vector<Node>& nodes, std::vector<u32>& mst, std::vector<Seq>& seqs,
                  u32* carry_out, u32* reps_out) {
  const u32 SEG = 1u << prm.lgseg;
  const u32 nice = prm.nice; const int lg = P.lgchunk;
  const int mode = P.litmode; const bool sub = lit_sub(mode); const int nctx = LIT_NCTX[mode];
  const u32 t2 = prm.far2, t3 = prm.far3;
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
      const u8* lmap = P.litmap + (size_t)P.ltype[ch] * nctx;
      const float* pl = P.ll[P.llmap[P.ctype[ch]]];
      const float* pm = P.mlk[P.mlkmap[P.ctype[ch] * MLK_NCTX + mlk_ctx(N.ll)]];
      const u8* omap = P.offmap + P.dtype[ch] * 4;
      const float* pal = P.align[P.almap[P.dtype[ch]]];
      u32 a = bs + sb + ss + i; const u8* cur = base + a; u32 maxl = std::min<u32>(sn - i, MAX_MATCH);
      u32 rl[3] = {0, 0, 0}; u32 bestr = 0, bestk = 0;
      if (maxl >= 2) for (u32 k = 0; k < 3; k++) {
        if ((k >= 1 && N.rep[k] == N.rep[0]) || (k == 2 && N.rep[2] == N.rep[1])) continue;
        rl[k] = match_len(cur, cur - N.rep[k], maxl);
        if (rl[k] > bestr) { bestr = rl[k]; bestk = k; }
      }
      const u32 m0 = mst[i], nm = mst[i + 1] - m0;
      const u32* mo = ms.off.data() + m0; const u16* mln = ms.len.data() + m0;
      u32 bestm = nm ? std::min<u32>(mln[nm - 1], maxl) : 0;
      const float base_c = N.cost + pl[0];
      auto match_price = [&](u32 ml, u32 kind) { u32 s, nb, ex; ml_code(ml - MIN_MATCH, s, nb, ex); return pm[kind * MLC_SYMS + s] + nb; };
      auto off_price = [&](u32 off, u32 c) {
        u32 s, nb, ex; off_code(off - 1, s, nb, ex);
        return P.off[omap[c]][s] + (nb >= 4 ? (float)(nb - 4) + pal[ex & 15] : (float)nb);
      };
      if (bestr >= nice || bestm >= nice) {
        if (bestr >= bestm) relax(i + bestr, base_c + match_price(bestr, bestk), N, bestr, N.rep[bestk], bestk);
        else relax(i + bestm, base_c + match_price(bestm, 3) + off_price(mo[nm - 1], 3), N, bestm, mo[nm - 1], 3);
        i += std::max(bestr, bestm) - 1;
        continue;
      }
      {  // literal
        u8 mb = *(cur - N.rep[0]);
        u32 ctx = lit_ctx(mode, a, cur[-1], cur[-2], mb);
        u8 s = sub ? (u8)(cur[0] - mb) : cur[0];
        float c = N.cost + P.lit[lmap[ctx]][s] + pll(pl, N.ll + 1) - pll(pl, N.ll);
        Node& T = nodes[i + 1];
        if (c < T.cost) { T = N; T.cost = c; T.ml = 0; T.ll = N.ll + 1; }
      }
      if (maxl < 2) continue;
      for (u32 k = 0; k < 3; k++) {
        for (u32 l = MIN_MATCH; l <= rl[k]; l++) relax(i + l, base_c + match_price(l, k), N, l, N.rep[k], k);
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
        for (u32 l = MIN_MATCH; l <= len; l++)
          relax(i + l, base_c + match_price(l, 3) + P.off[omap[std::min<u32>(l - MIN_MATCH, 3)]][OFF_SLOTS + sp], N, l, off, 4 + sp);
      }
      u32 prev = MIN_MATCH - 1;
      for (u32 j = 0; j < nm; j++) {
        u32 lc = std::min<u32>(mln[j], maxl); u32 off = mo[j];
        if (lc <= prev) continue;
        if (off == N.rep[0] || off == N.rep[1] || off == N.rep[2]) { prev = lc; continue; }
        u32 lmin = off > t3 ? 4 : off > t2 ? 3 : 2;
        float op[4]; for (u32 c = 0; c < 4; c++) op[c] = off_price(off, c);
        for (u32 l = std::max(prev + 1, lmin); l <= lc; l++)
          relax(i + l, base_c + match_price(l, 3) + op[std::min<u32>(l - MIN_MATCH, 3)], N, l, off, 3);
        prev = lc;
      }
    }
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
struct Ev { u16 start, freq; u8 st; };
static double g_bits[8];  // lit, ll, mlk, off, align, len-extra, off-extra
static u64 g_kinds[20];
static u64 g_tables[8];

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
  auto add = [&](u16 start, u16 freq, u8 st) {
    if (!freq) { fprintf(stderr, "internal error: zero-frequency symbol\n"); abort(); }
    ev.push_back(Ev{start, freq, st});
  };
  u32 pos = bs, li = 0; u32 r[NUM_REPS]; memcpy(r, reps_in, sizeof r);
  const int mode = M.litmode; const bool sub = lit_sub(mode); const int nctx = LIT_NCTX[mode]; const int lg = M.lgchunk;
  const u32 end = bs + bn; u32 chunk_end = bs;
  const u8* lmap = nullptr; const u8* mmap = nullptr; const u8* omap = nullptr; u32 llc = 0, alc = 0;
  auto upd = [&]() {
    if (pos >= end) return;
    u32 ch = (pos - bs) >> lg; chunk_end = bs + ((ch + 1) << lg);
    lmap = M.litmap + (size_t)M.ltype[ch] * nctx; llc = M.llmap[M.ctype[ch]];
    mmap = M.mlkmap + M.ctype[ch] * MLK_NCTX; omap = M.offmap + M.dtype[ch] * 4; alc = M.almap[M.dtype[ch]];
  };
  auto lits = [&](u32 cnt) {
    while (cnt) {
      u32 n = std::min(cnt, chunk_end - pos);
      for (u32 k = 0; k < n; k++, pos++) {
        const u8* q = base + pos; u8 mb = *(q - r[0]); u32 cl = lmap[lit_ctx(mode, pos, *(q - 1), *(q - 2), mb)];
        u8 s = sub ? (u8)(base[pos] - mb) : base[pos];
        add(clit[cl][s], M.lit[cl][s], (u8)(li & 1)); li++; g_bits[0] += FPRICE[M.lit[cl][s]];
      }
      cnt -= n; if (pos >= chunk_end) upd();
    }
  };
  for (const Seq& s : seqs) {
    u32 sym, nb, ex;
    if (pos >= chunk_end) upd();
    ll_code(s.ll, sym, nb, ex); add(cll[llc][sym], M.ll[llc][sym], 2); bits.put(ex, nb); g_bits[1] += FPRICE[M.ll[llc][sym]]; g_bits[5] += nb;
    lits(s.ll);
    ml_code(s.ml - MIN_MATCH, sym, nb, ex);
    u32 c = mmap[mlk_ctx(s.ll)]; u32 ms = std::min<u32>(s.kind, 3) * MLC_SYMS + sym;
    add(cmlk[c][ms], M.mlk[c][ms], 3); bits.put(ex, nb); g_bits[2] += FPRICE[M.mlk[c][ms]]; g_bits[5] += nb;
    g_kinds[s.kind]++;
    if (s.kind >= 3) {
      sym = off_symbol(s, nb, ex);
      u32 oc = omap[std::min<u32>(s.ml - MIN_MATCH, 3)];
      add(coff[oc][sym], M.off[oc][sym], 2); g_bits[3] += FPRICE[M.off[oc][sym]]; g_bits[6] += nb >= 4 ? nb - 4 : nb;
      if (nb >= 4) { bits.put(ex >> 4, nb - 4); add(cal[alc][ex & 15], M.align[alc][ex & 15], 3); g_bits[4] += FPRICE[M.align[alc][ex & 15]]; }
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
    const Ev& e = ev[i]; u32 xs = x[e.st];
    u64 xmax = (u64)((RANS_L >> SCALE_BITS) << 16) * e.freq;
    if (xs >= xmax) { *--ptr = (u16)xs; xs >>= 16; }
    x[e.st] = ((xs / e.freq) << SCALE_BITS) + (xs % e.freq) + e.start;
  }
  for (int i = NUM_STREAMS - 1; i >= 0; i--) { *--ptr = (u16)(x[i] >> 16); *--ptr = (u16)x[i]; }
  size_t rb = (w.data() + w.size() - ptr) * 2;
  // header
  BitWriter hw;
  hw.put(lg, 5);
  hw.put(M.nlt - 1, 4); write_map(hw, M.ltype, M.nch, M.nlt);
  hw.put(M.nct - 1, 4); write_map(hw, M.ctype, M.nch, M.nct);
  hw.put(M.ndt - 1, 4); write_map(hw, M.dtype, M.nch, M.ndt);
  hw.put(mode, 4); hw.put(M.ncl - 1, 6); write_map(hw, M.litmap, M.nlt * nctx, M.ncl);
  for (int c = 0; c < M.ncl; c++) write_table(hw, M.lit[c], 256);
  hw.put(M.nll - 1, 3); write_map(hw, M.llmap, M.nct, M.nll);
  for (int c = 0; c < M.nll; c++) write_table(hw, M.ll[c], LL_SYMS);
  hw.put(M.nmk - 1, 4); write_map(hw, M.mlkmap, M.nct * MLK_NCTX, M.nmk);
  for (int c = 0; c < M.nmk; c++) write_table(hw, M.mlk[c], MLK_SYMS);
  hw.put(M.nof - 1, 4); write_map(hw, M.offmap, M.ndt * 4, M.nof);
  for (int c = 0; c < M.nof; c++) write_table(hw, M.off[c], OFF_SYMS);
  hw.put(M.nal - 1, 2); write_map(hw, M.almap, M.ndt, M.nal);
  for (int c = 0; c < M.nal; c++) write_table(hw, M.align[c], ALIGN_SYMS);
  hw.flush();
  g_tables[0] += M.ncl; g_tables[1] += M.nll; g_tables[2] += M.nmk; g_tables[3] += M.nof; g_tables[4] += M.nlt; g_tables[5] += M.nct; g_tables[6] += M.ndt; g_tables[7]++;
  put_varint(out, bn); out.push_back(1);
  put_varint(out, seqs.size()); put_varint(out, hw.buf.size()); put_varint(out, rb); put_varint(out, bits.buf.size());
  out.insert(out.end(), hw.buf.begin(), hw.buf.end());
  out.insert(out.end(), (const u8*)ptr, (const u8*)ptr + rb);
  out.insert(out.end(), bits.buf.begin(), bits.buf.end());
  if (hdr_bytes) *hdr_bytes += hw.buf.size();
}

std::vector<u8> compress(const u8* base, size_t n, const Params& prm) {
  einit();
  std::vector<u8> out = {'G', 'T', '0', 2};
  put_varint(out, n);
  MF mf(base, (u32)n, prm);
  Prices* pr = new Prices; default_prices(*pr, prm.lgchunk);
  Models* M = new Models;
  u32 reps[NUM_REPS] = {1, 2, 3, 4};
  const u32 B = 1u << prm.lgblock;
  const u32 SB = std::min(B, 1u << prm.lgsub);
  MatchStore ms; std::vector<u32> mst; std::vector<Match> tmp(prm.nice + 8);
  std::vector<Node> nodes; std::vector<Seq> seqs, sseqs; LitInfo L;
  int modestat[NUM_LITMODES] = {0}; size_t hdr_bytes = 0, nseq_tot = 0, stored = 0;
  bool prior = false;
  const u32 t2 = prm.far2, t3 = prm.far3;
  for (u64 bs64 = 0; bs64 < n; bs64 += B) {
    const u32 bs = (u32)bs64; const u32 bn = (u32)std::min<u64>(B, n - bs64);
    carry_prior(*pr, (int)((bn + (1u << prm.lgchunk) - 1) >> prm.lgchunk));
    seqs.clear();
    u32 sreps[NUM_REPS]; memcpy(sreps, reps, sizeof sreps); u32 carry = 0;
    for (u32 sb = 0; sb < bn; sb += SB) {
      const u32 sn = std::min(SB, bn - sb);
      double t0 = now();
      ms.cnt.assign(sn, 0); ms.off.clear(); ms.len.clear();
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
      int passes = prior ? prm.passes : prm.passes + 1;
      u32 creps[NUM_REPS], cout = 0;
      for (int it = 0; it < passes; it++) {
        t0 = now();
        parse(base, bs, sb, sn, sreps, carry, ms, *pr, prm, nodes, mst, sseqs, &cout, creps);
        g_time[1] += now() - t0; t0 = now();
        if (SB < bn || it < passes - 1) {
          build_models(base, bs, bn, bs + sb - carry, bs + sb + sn, sreps, sseqs, prm, *M, L, nullptr, prm.fullevery || (sb == 0 && it == 0));
          prices_from_models(*M, *pr);
        }
        g_time[2] += now() - t0;
      }
      prior = true;
      seqs.insert(seqs.end(), sseqs.begin(), sseqs.end());
      carry = cout; memcpy(sreps, creps, sizeof sreps);
    }
    double t1 = now();
    build_models(base, bs, bn, bs, bs + bn, reps, seqs, prm, *M, L, modestat, true);
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
  }
  if (prm.verbose) {
    fprintf(stderr, "time: matchfind %.2f parse %.2f models %.2f emit %.2f\n", g_time[0], g_time[1], g_time[2], g_time[3]);
    fprintf(stderr, "litmodes:"); for (int i = 0; i < NUM_LITMODES; i++) fprintf(stderr, " %d", modestat[i]);
    fprintf(stderr, "  stored %zu  hdr %zu B  seqs %zu\n", stored, hdr_bytes, nseq_tot);
    double nb = (double)std::max<u64>(1, g_tables[7]);
    fprintf(stderr, "avg/block: lit tables %.1f ll %.1f mlk %.1f off %.1f | types lit %.1f cmd %.1f dist %.1f\n", g_tables[0] / nb, g_tables[1] / nb, g_tables[2] / nb, g_tables[3] / nb, g_tables[4] / nb, g_tables[5] / nb, g_tables[6] / nb);
    fprintf(stderr, "bytes: lit %.0f ll %.0f mlk %.0f off %.0f align %.0f len-extra %.0f off-extra %.0f\n", g_bits[0] / 8, g_bits[1] / 8, g_bits[2] / 8, g_bits[3] / 8, g_bits[4] / 8, g_bits[5] / 8, g_bits[6] / 8);
    fprintf(stderr, "kinds: rep0 %llu rep1 %llu rep2 %llu new %llu rep3 %llu near0", (unsigned long long)g_kinds[0], (unsigned long long)g_kinds[1], (unsigned long long)g_kinds[2], (unsigned long long)g_kinds[3], (unsigned long long)g_kinds[4]);
    u64 a = 0, b = 0; for (int i = 5; i < 11; i++) a += g_kinds[i]; for (int i = 11; i < 17; i++) b += g_kinds[i];
    fprintf(stderr, " %llu near1 %llu\n", (unsigned long long)a, (unsigned long long)b);
  }
  delete pr; delete M;
  return out;
}

}  // namespace gt
