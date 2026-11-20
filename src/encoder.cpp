#include "common.h"
namespace gt {

static const u32 NIL = 0xffffffffu;
static const float INF = 1e30f;

struct Match { u32 len, off; };
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
    head.assign((size_t)1 << hb, NIL); son.resize(2 * (size_t)W);
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
      out[cnt - 1].len = match_len(cur, cur - out[cnt - 1].off, std::min<u32>(avail, MAX_MATCH));
    return cnt;
  }
};

// ---------------- models / prices ----------------
struct Models {
  int litmode, ncl; u8 map[MAXCTX];
  u16 lit[MAXCL][256]; u16 ll[LL_SYMS]; u16 mlk[MLK_MAXCL][MLK_SYMS]; u16 off[4][OFF_SYMS]; u16 align[ALIGN_SYMS];
  int nmk; u8 mlkmap[MLK_NCTX];
};
struct Prices {
  int litmode; u8 map[MAXCTX];
  float lit[MAXCL][256]; float ll[LL_SYMS]; float mlk[MLK_MAXCL][MLK_SYMS]; float off[4][OFF_SYMS]; float align[ALIGN_SYMS];
  u8 mlkmap[MLK_NCTX];
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

static void default_prices(Prices& P) {
  P.litmode = 0; memset(P.map, 0, sizeof P.map);
  for (int c = 0; c < MAXCL; c++) for (int s = 0; s < 256; s++) P.lit[c][s] = 6.5f;
  for (int s = 0; s < LL_SYMS; s++) P.ll[s] = 2.0f + 0.25f * s;
  memset(P.mlkmap, 0, sizeof P.mlkmap);
  for (int c = 0; c < 1; c++) for (int s = 0; s < MLK_SYMS; s++) {
    int k = s / MLC_SYMS, mc = s % MLC_SYMS;
    P.mlk[c][s] = (k == 0 ? 3.0f : k == 3 ? 4.0f : 6.0f) + 0.1f * std::min(mc, 32);
  }
  for (int c = 0; c < 4; c++) for (int s = 0; s < OFF_SYMS; s++) P.off[c][s] = 5.5f;
  for (int s = 0; s < ALIGN_SYMS; s++) P.align[s] = 4.0f;
}
static void prices_from_models(const Models& M, Prices& P) {
  P.litmode = M.litmode; memcpy(P.map, M.map, sizeof P.map);
  for (int c = 0; c < M.ncl; c++) for (int s = 0; s < 256; s++) P.lit[c][s] = FPRICE[M.lit[c][s]];
  for (int s = 0; s < LL_SYMS; s++) P.ll[s] = FPRICE[M.ll[s]];
  memcpy(P.mlkmap, M.mlkmap, sizeof P.mlkmap);
  for (int c = 0; c < M.nmk; c++) for (int s = 0; s < MLK_SYMS; s++) P.mlk[c][s] = FPRICE[M.mlk[c][s]];
  for (int c = 0; c < 4; c++) for (int s = 0; s < OFF_SYMS; s++) P.off[c][s] = FPRICE[M.off[c][s]];
  for (int s = 0; s < ALIGN_SYMS; s++) P.align[s] = FPRICE[M.align[s]];
}

// ---------------- context clustering ----------------
static inline double hdr_est(int nz) { return 32 + 8.0 * nz; }
static double hist_cost(const Hist& h) {
  double s = xlogx(h.tot); int nz = 0;
  for (int i = 0; i < 256; i++) if (h.c[i]) { s -= xlogx(h.c[i]); nz++; }
  return s + hdr_est(nz);
}
// Greedy agglomerative merge of cl[first..first+n) (all alive). Merges while it reduces cost or
// while more than maxcl remain. merged_into[i] gets the surviving index for every input.
static void merge_range(std::vector<Hist>& cl, std::vector<double>& cc, int first, int n, int maxcl, std::vector<int>& into) {
  std::vector<double> D((size_t)n * n, 0.0);
  std::vector<char> alive(n, 1);
  auto pd = [&](int a, int b) {
    const Hist& A = cl[first + a]; const Hist& B = cl[first + b];
    double s = xlogx(A.tot + B.tot); int nz = 0;
    for (int i = 0; i < 256; i++) { u32 v = A.c[i] + B.c[i]; if (v) { s -= xlogx(v); nz++; } }
    return s + hdr_est(nz) - cc[first + a] - cc[first + b];
  };
  for (int a = 0; a < n; a++) { into[first + a] = first + a; for (int b = a + 1; b < n; b++) D[(size_t)a * n + b] = pd(a, b); }
  int nalive = n;
  while (nalive > 1) {
    double best = 1e300; int ba = -1, bb = -1;
    for (int a = 0; a < n; a++) if (alive[a]) for (int b = a + 1; b < n; b++) if (alive[b]) {
      double d = D[(size_t)a * n + b]; if (d < best) { best = d; ba = a; bb = b; }
    }
    if (best >= 0 && nalive <= maxcl) break;
    Hist& A = cl[first + ba]; const Hist& B = cl[first + bb];
    for (int i = 0; i < 256; i++) A.c[i] += B.c[i];
    A.tot += B.tot; cc[first + ba] = hist_cost(A); alive[bb] = 0; nalive--;
    for (int j = 0; j < n; j++) if (into[first + j] == first + bb) into[first + j] = first + ba;
    for (int j = 0; j < n; j++) if (alive[j] && j != ba) {
      int a = std::min(j, ba), b = std::max(j, ba); D[(size_t)a * n + b] = pd(a, b);
    }
  }
}
// returns number of clusters; fills map[ctx] and out[]; cost = estimated total bits
static int cluster(const std::vector<Hist>& h, int nctx, int maxcl, u8* map, std::vector<Hist>& out, double& cost) {
  std::vector<int> owner(nctx, -1);
  std::vector<Hist> cl; std::vector<double> cc;
  for (int i = 0; i < nctx; i++) if (h[i].tot) { owner[i] = (int)cl.size(); cl.push_back(h[i]); cc.push_back(hist_cost(h[i])); }
  int n = (int)cl.size();
  out.clear();
  if (n == 0) { memset(map, 0, nctx); Hist z; memset(&z, 0, sizeof z); out.push_back(z); cost = 0; return 1; }
  std::vector<int> into(n);
  if (n > 96) {  // stage 1: reduce in chunks, then compact survivors
    for (int f = 0; f < n; f += 64) merge_range(cl, cc, f, std::min(64, n - f), 12, into);
    std::vector<int> newidx(n, -1); std::vector<Hist> c2; std::vector<double> cc2;
    for (int i = 0; i < n; i++) if (into[i] == i) { newidx[i] = (int)c2.size(); c2.push_back(cl[i]); cc2.push_back(cc[i]); }
    for (int i = 0; i < nctx; i++) if (owner[i] >= 0) owner[i] = newidx[into[owner[i]]];
    cl.swap(c2); cc.swap(cc2); n = (int)cl.size(); into.assign(n, 0);
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
struct LitInfo { std::vector<u8> v, p1, p2, mb, pos; };

static inline void rep_update(u32* r, u32 kind, u32 off) {
  if (kind == 0) return;
  if (kind == 1) { r[1] = r[0]; r[0] = off; return; }
  if (kind == 2) { r[2] = r[1]; r[1] = r[0]; r[0] = off; return; }
  r[3] = r[2]; r[2] = r[1]; r[1] = r[0]; r[0] = off;
}

static void replay(const u8* base, u32 bs, u32 bn, const u32* reps_in, const std::vector<Seq>& seqs, LitInfo* L, u32* reps_out) {
  u32 pos = bs; u32 r[NUM_REPS]; memcpy(r, reps_in, sizeof r);
  if (L) { L->v.clear(); L->p1.clear(); L->p2.clear(); L->mb.clear(); L->pos.clear(); }
  auto lits = [&](u32 cnt) {
    if (L) for (u32 k = 0; k < cnt; k++, pos++) {
      const u8* q = base + pos;
      L->v.push_back(q[0]); L->p1.push_back(*(q - 1)); L->p2.push_back(*(q - 2)); L->mb.push_back(*(q - r[0])); L->pos.push_back((u8)(pos & 3));
    } else pos += cnt;
  };
  for (const Seq& s : seqs) { lits(s.ll); rep_update(r, s.kind, s.off); pos += s.ml; }
  lits(bs + bn - pos);
  if (reps_out) memcpy(reps_out, r, sizeof r);
}

static inline u32 off_symbol(const Seq& s, u32& nb, u32& ex) {  // for kind >= 3
  if (s.kind == 3) { u32 sym; off_code(s.off - 1, sym, nb, ex); return sym; }
  nb = 0; ex = 0; return OFF_SLOTS + (s.kind - 4);
}

static void build_models(const u8* base, u32 bs, u32 bn, const u32* reps_in, const std::vector<Seq>& seqs,
                         const Params& prm, Models& M, LitInfo& L, int* modestat) {
  replay(base, bs, bn, reps_in, seqs, &L, nullptr);
  size_t nl = L.v.size();
  double bestc = 1e300; std::vector<Hist> besth; static u8 bestmap[MAXCTX], map[MAXCTX]; int bestmode = 0;
  std::vector<Hist> h, outh;
  for (int mode = 0; mode < NUM_LITMODES; mode++) {
    if (!((prm.modemask >> mode) & 1)) continue;
    int nctx = LIT_NCTX[mode];
    h.resize(nctx); memset(h.data(), 0, sizeof(Hist) * nctx);
    bool sub = lit_sub(mode);
    for (size_t i = 0; i < nl; i++) {
      u32 ctx = lit_ctx(mode, L.pos[i], L.p1[i], L.p2[i], L.mb[i]);
      u8 s = sub ? (u8)(L.v[i] - L.mb[i]) : L.v[i];
      h[ctx].c[s]++; h[ctx].tot++;
    }
    double cost; cluster(h, nctx, prm.maxcl, map, outh, cost);
    if (cost < bestc) { bestc = cost; besth = outh; memcpy(bestmap, map, sizeof map); bestmode = mode; }
  }
  M.litmode = bestmode; M.ncl = (int)besth.size(); memcpy(M.map, bestmap, sizeof bestmap);
  for (int c = 0; c < M.ncl; c++) normalize_freqs(besth[c].c, 256, M.lit[c]);
  if (modestat) modestat[bestmode]++;
  u32 hll[LL_SYMS] = {0}, hoff[4][OFF_SYMS] = {{0}}, hal[ALIGN_SYMS] = {0};
  std::vector<Hist> hmlk(MLK_NCTX); memset(hmlk.data(), 0, sizeof(Hist) * MLK_NCTX);
  for (const Seq& s : seqs) {
    u32 sym, nb, ex;
    ll_code(s.ll, sym, nb, ex); hll[sym]++;
    ml_code(s.ml - MIN_MATCH, sym, nb, ex); { Hist& hh = hmlk[mlk_ctx(s.ll)]; hh.c[std::min<u32>(s.kind, 3) * MLC_SYMS + sym]++; hh.tot++; }
    if (s.kind >= 3) {
      sym = off_symbol(s, nb, ex); hoff[std::min<u32>(s.ml - MIN_MATCH, 3)][sym]++;
      if (nb >= 4) hal[ex & 15]++;
    }
  }
  normalize_freqs(hll, LL_SYMS, M.ll);
  { double cost; std::vector<Hist> o; M.nmk = cluster(hmlk, MLK_NCTX, MLK_MAXCL, M.mlkmap, o, cost);
    for (int c = 0; c < M.nmk; c++) normalize_freqs(o[c].c, MLK_SYMS, M.mlk[c]); }
  for (int c = 0; c < 4; c++) normalize_freqs(hoff[c], OFF_SYMS, M.off[c]);
  normalize_freqs(hal, ALIGN_SYMS, M.align);
}

// ---------------- optimal parse ----------------
struct Node { float cost; u32 ml, off, ll, kind; u32 rep[NUM_REPS]; };

static inline float pll(const Prices& P, u32 v) { u32 s, nb, ex; ll_code(v, s, nb, ex); return P.ll[s] + nb; }

static void parse(const u8* base, u32 bs, u32 bn, const u32* reps_in, const std::vector<u32>& mstart,
                  const std::vector<Match>& mlist, const Prices& P, const Params& prm,
                  std::vector<Node>& nodes, std::vector<Seq>& seqs) {
  nodes.resize(bn + 1);
  for (u32 i = 0; i <= bn; i++) nodes[i].cost = INF;
  Node& n0 = nodes[0]; n0.cost = 0; n0.ml = 0; n0.off = 0; n0.ll = 0; n0.kind = 0;
  memcpy(n0.rep, reps_in, sizeof n0.rep);
  const float ll0 = pll(P, 0);
  const u32 nice = prm.nice;
  const int mode = P.litmode; const bool sub = lit_sub(mode);
  auto relax = [&](u32 t, float c, const Node& N, u32 ml, u32 off, u32 kind) {
    Node& T = nodes[t];
    if (c < T.cost) {
      T.cost = c; T.ml = ml; T.off = off; T.kind = kind; T.ll = 0;
      memcpy(T.rep, N.rep, sizeof T.rep); rep_update(T.rep, kind, off);
    }
  };
  for (u32 i = 0; i < bn; i++) {
    const Node N = nodes[i];
    if (N.cost >= INF) continue;
    u32 a = bs + i; const u8* cur = base + a; u32 maxl = std::min<u32>(bn - i, MAX_MATCH);
    u32 rl[3] = {0, 0, 0}; u32 bestr = 0, bestk = 0;
    if (maxl >= 2) for (u32 k = 0; k < 3; k++) {
      if ((k >= 1 && N.rep[k] == N.rep[0]) || (k == 2 && N.rep[2] == N.rep[1])) continue;
      rl[k] = match_len(cur, cur - N.rep[k], maxl);
      if (rl[k] > bestr) { bestr = rl[k]; bestk = k; }
    }
    const Match* m = mlist.data() + mstart[i]; u32 nm = mstart[i + 1] - mstart[i];
    u32 bestm = nm ? std::min(m[nm - 1].len, maxl) : 0;
    const float* pm = P.mlk[P.mlkmap[mlk_ctx(N.ll)]];
    const float base_c = N.cost + ll0;
    auto match_price = [&](u32 ml, u32 kind) { u32 s, nb, ex; ml_code(ml - MIN_MATCH, s, nb, ex); return pm[kind * MLC_SYMS + s] + nb; };
    auto off_price = [&](u32 off, u32 c) {
      u32 s, nb, ex; off_code(off - 1, s, nb, ex);
      return P.off[c][s] + (nb >= 4 ? (float)(nb - 4) + P.align[ex & 15] : (float)nb);
    };
    if (bestr >= nice || bestm >= nice) {
      if (bestr >= bestm) relax(i + bestr, base_c + match_price(bestr, bestk), N, bestr, N.rep[bestk], bestk);
      else relax(i + bestm, base_c + match_price(bestm, 3) + off_price(m[nm - 1].off, 3), N, bestm, m[nm - 1].off, 3);
      i += std::max(bestr, bestm) - 1;
      continue;
    }
    {  // literal
      u8 mb = *(cur - N.rep[0]);
      u32 ctx = lit_ctx(mode, a, cur[-1], cur[-2], mb);
      u8 s = sub ? (u8)(cur[0] - mb) : cur[0];
      float c = N.cost + P.lit[P.map[ctx]][s] + pll(P, N.ll + 1) - pll(P, N.ll);
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
        relax(i + l, base_c + match_price(l, 3) + P.off[std::min<u32>(l - MIN_MATCH, 3)][OFF_SLOTS + sp], N, l, off, 4 + sp);
    }
    u32 prev = MIN_MATCH - 1;
    for (u32 j = 0; j < nm; j++) {
      u32 lc = std::min(m[j].len, maxl); u32 off = m[j].off;
      if (lc <= prev) continue;
      if (off == N.rep[0] || off == N.rep[1] || off == N.rep[2]) { prev = lc; continue; }
      float op[4]; for (u32 c = 0; c < 4; c++) op[c] = off_price(off, c);
      for (u32 l = prev + 1; l <= lc; l++)
        relax(i + l, base_c + match_price(l, 3) + op[std::min<u32>(l - MIN_MATCH, 3)], N, l, off, 3);
      prev = lc;
    }
  }
  std::vector<Seq> rev; u32 run = 0; u32 i = bn;
  while (i > 0) {
    const Node& n = nodes[i];
    if (n.ml == 0) { run++; i--; }
    else {
      if (!rev.empty()) rev.back().ll = run;
      rev.push_back(Seq{0, n.ml, n.off, n.kind}); run = 0; i -= n.ml;
    }
  }
  if (!rev.empty()) rev.back().ll = run;
  seqs.assign(rev.rbegin(), rev.rend());
}

// ---------------- block emission ----------------
struct Ev { u16 start, freq; u8 st; };
static double g_bits[8];  // lit, ll, mlk, off, align, extra
static u64 g_kinds[20];

static void cumul(const u16* f, int n, u16* cum) { u32 c = 0; for (int i = 0; i < n; i++) { cum[i] = (u16)c; c += f[i]; } }

static void emit_block(const u8* base, u32 bs, u32 bn, const u32* reps_in, const std::vector<Seq>& seqs,
                       const Models& M, std::vector<u8>& out, size_t* hdr_bytes) {
  static u16 clit[MAXCL][256]; u16 cll[LL_SYMS], cmlk[MLK_MAXCL][MLK_SYMS], coff[4][OFF_SYMS], cal[ALIGN_SYMS];
  for (int c = 0; c < M.ncl; c++) cumul(M.lit[c], 256, clit[c]);
  cumul(M.ll, LL_SYMS, cll);
  for (int c = 0; c < M.nmk; c++) cumul(M.mlk[c], MLK_SYMS, cmlk[c]);
  for (int c = 0; c < 4; c++) cumul(M.off[c], OFF_SYMS, coff[c]);
  cumul(M.align, ALIGN_SYMS, cal);
  std::vector<Ev> ev; ev.reserve(bn + 4 * seqs.size() + 16);
  BitWriter bits;
  auto add = [&](u16 start, u16 freq, u8 st) {
    if (!freq) { fprintf(stderr, "internal error: zero-frequency symbol\n"); abort(); }
    ev.push_back(Ev{start, freq, st});
  };
  u32 pos = bs, li = 0; u32 r[NUM_REPS]; memcpy(r, reps_in, sizeof r);
  const int mode = M.litmode; const bool sub = lit_sub(mode);
  auto lits = [&](u32 cnt) {
    for (u32 k = 0; k < cnt; k++, pos++) {
      const u8* q = base + pos; u8 mb = *(q - r[0]); u32 cl = M.map[lit_ctx(mode, pos, *(q - 1), *(q - 2), mb)];
      u8 s = sub ? (u8)(base[pos] - mb) : base[pos];
      add(clit[cl][s], M.lit[cl][s], (u8)(li & 1)); li++; g_bits[0] += FPRICE[M.lit[cl][s]];
    }
  };
  for (const Seq& s : seqs) {
    u32 sym, nb, ex;
    ll_code(s.ll, sym, nb, ex); add(cll[sym], M.ll[sym], 2); bits.put(ex, nb); g_bits[1] += FPRICE[M.ll[sym]]; g_bits[5] += nb;
    lits(s.ll);
    ml_code(s.ml - MIN_MATCH, sym, nb, ex);
    int c = M.mlkmap[mlk_ctx(s.ll)]; u32 ms = std::min<u32>(s.kind, 3) * MLC_SYMS + sym;
    add(cmlk[c][ms], M.mlk[c][ms], 3); bits.put(ex, nb); g_bits[2] += FPRICE[M.mlk[c][ms]]; g_bits[5] += nb;
    g_kinds[s.kind]++;
    if (s.kind >= 3) {
      sym = off_symbol(s, nb, ex);
      u32 oc = std::min<u32>(s.ml - MIN_MATCH, 3);
      add(coff[oc][sym], M.off[oc][sym], 2); g_bits[3] += FPRICE[M.off[oc][sym]]; g_bits[6] += nb >= 4 ? nb - 4 : nb;
      if (nb >= 4) { bits.put(ex >> 4, nb - 4); add(cal[ex & 15], M.align[ex & 15], 3); g_bits[4] += FPRICE[M.align[ex & 15]]; }
      else bits.put(ex, nb);
    }
    rep_update(r, s.kind, s.off);
    pos += s.ml;
  }
  lits(bs + bn - pos);
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
  hw.put(mode, 4); hw.put(M.ncl, 7);
  write_map(hw, M.map, LIT_NCTX[mode], M.ncl);
  for (int c = 0; c < M.ncl; c++) write_table(hw, M.lit[c], 256);
  write_table(hw, M.ll, LL_SYMS);
  hw.put(M.nmk - 1, 3); write_map(hw, M.mlkmap, MLK_NCTX, M.nmk);
  for (int c = 0; c < M.nmk; c++) write_table(hw, M.mlk[c], MLK_SYMS);
  for (int c = 0; c < 4; c++) write_table(hw, M.off[c], OFF_SYMS);
  write_table(hw, M.align, ALIGN_SYMS);
  hw.flush();
  put_varint(out, bn); out.push_back(1);
  put_varint(out, seqs.size()); put_varint(out, hw.buf.size()); put_varint(out, rb); put_varint(out, bits.buf.size());
  out.insert(out.end(), hw.buf.begin(), hw.buf.end());
  out.insert(out.end(), (const u8*)ptr, (const u8*)ptr + rb);
  out.insert(out.end(), bits.buf.begin(), bits.buf.end());
  if (hdr_bytes) *hdr_bytes += hw.buf.size();
}

std::vector<u8> compress(const u8* in, size_t n, const Params& prm) {
  einit();
  std::vector<u8> buf(n + 2 * PAD + 64, 0);
  if (n) memcpy(buf.data() + PAD, in, n);
  const u8* base = buf.data() + PAD;
  std::vector<u8> out = {'G', 'T', '0', 1};
  put_varint(out, n);
  MF mf(base, (u32)n, prm);
  Prices* pr = new Prices; default_prices(*pr);
  Models* M = new Models;
  u32 reps[NUM_REPS] = {1, 2, 3, 4};
  const u32 B = 1u << prm.lgblock;
  std::vector<u32> mstart; std::vector<Match> mlist; std::vector<Match> tmp(prm.nice + 8);
  std::vector<Node> nodes; std::vector<Seq> seqs; LitInfo L;
  int modestat[NUM_LITMODES] = {0}; size_t hdr_bytes = 0, nseq_tot = 0, stored = 0;
  bool prior = false;
  for (u64 bs = 0; bs < n; bs += B) {
    u32 bn = (u32)std::min<u64>(B, n - bs);
    mstart.assign(bn + 1, 0); mlist.clear();
    for (u32 i = 0; i < bn; i++) {
      mstart[i] = (u32)mlist.size();
      u32 cnt = mf.find((u32)bs + i, tmp.data(), true);
      mlist.insert(mlist.end(), tmp.begin(), tmp.begin() + cnt);
      if (cnt && tmp[cnt - 1].len >= (u32)prm.nice) {
        u32 Lm = std::min(tmp[cnt - 1].len, bn - i);
        for (u32 j = 1; j < Lm; j++) { mstart[i + j] = (u32)mlist.size(); mf.find((u32)bs + i + j, tmp.data(), false); }
        i += Lm - 1;
      }
    }
    mstart[bn] = (u32)mlist.size();
    int passes = prior ? prm.passes : prm.passes + 1;
    for (int it = 0; it < passes; it++) {
      parse(base, (u32)bs, bn, reps, mstart, mlist, *pr, prm, nodes, seqs);
      build_models(base, (u32)bs, bn, reps, seqs, prm, *M, L, it == passes - 1 ? modestat : nullptr);
      prices_from_models(*M, *pr);
    }
    prior = true;
    size_t before = out.size();
    emit_block(base, (u32)bs, bn, reps, seqs, *M, out, &hdr_bytes);
    if (out.size() - before >= bn + 6) {
      out.resize(before); put_varint(out, bn); out.push_back(0);
      out.insert(out.end(), base + bs, base + bs + bn); stored++;
    } else {
      replay(base, (u32)bs, bn, reps, seqs, nullptr, reps);
      nseq_tot += seqs.size();
    }
  }
  if (prm.verbose) {
    fprintf(stderr, "litmodes:"); for (int i = 0; i < NUM_LITMODES; i++) fprintf(stderr, " %d", modestat[i]);
    fprintf(stderr, "  stored %zu  hdr %zu B  seqs %zu\n", stored, hdr_bytes, nseq_tot);
    fprintf(stderr, "bytes: lit %.0f ll %.0f mlk %.0f off %.0f align %.0f len-extra %.0f off-extra %.0f\n", g_bits[0] / 8, g_bits[1] / 8, g_bits[2] / 8, g_bits[3] / 8, g_bits[4] / 8, g_bits[5] / 8, g_bits[6] / 8);
    fprintf(stderr, "kinds: rep0 %llu rep1 %llu rep2 %llu new %llu rep3 %llu near0", (unsigned long long)g_kinds[0], (unsigned long long)g_kinds[1], (unsigned long long)g_kinds[2], (unsigned long long)g_kinds[3], (unsigned long long)g_kinds[4]);
    u64 a = 0, b = 0; for (int i = 5; i < 11; i++) a += g_kinds[i]; for (int i = 11; i < 17; i++) b += g_kinds[i];
    fprintf(stderr, " %llu near1 %llu\n", (unsigned long long)a, (unsigned long long)b);
  }
  delete pr; delete M;
  return out;
}

}  // namespace gt
