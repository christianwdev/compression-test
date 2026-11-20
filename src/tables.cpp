#include "common.h"
namespace gt {

void normalize_freqs(const u32* cnt, int n, u16* f, int bits) {
  const int SCALE = 1 << bits;
  u64 tot = 0;
  for (int i = 0; i < n; i++) tot += cnt[i];
  if (tot == 0) { memset(f, 0, n * sizeof(u16)); f[0] = SCALE - 1; f[1] = 1; return; }
  int sum = 0;
  for (int i = 0; i < n; i++) {
    if (!cnt[i]) { f[i] = 0; continue; }
    u64 v = ((u64)cnt[i] * SCALE + tot / 2) / tot;
    if (v < 1) v = 1;
    if (v > (u64)SCALE - 1) v = (u64)SCALE - 1;
    f[i] = (u16)v; sum += (int)v;
  }
  int diff = SCALE - sum;
  while (diff > 0) {
    int best = -1; double bg = -1;
    for (int i = 0; i < n; i++) if (cnt[i] && f[i] < SCALE - 1) {
      double g = cnt[i] * log2((f[i] + 1.0) / f[i]);
      if (g > bg) { bg = g; best = i; }
    }
    if (best < 0) { for (int i = 0; i < n; i++) if (!f[i]) { best = i; break; } }
    f[best]++; diff--;
  }
  while (diff < 0) {
    int best = -1; double bl = 1e300;
    for (int i = 0; i < n; i++) if (f[i] > 1) {
      double l = cnt[i] * log2(f[i] / (f[i] - 1.0));
      if (l < bl) { bl = l; best = i; }
    }
    f[best]--; diff++;
  }
}

// ---- adaptive binary range coder (the LZMA construction), used only for block headers
enum { RC_BITS = 11, RC_MOVE = 4, RC_TOP = 1 << 24 };

void RcEnc::shift_low() {
  if ((u32)low < 0xff000000u || (low >> 32) != 0) {
    u8 carry = (u8)(low >> 32), c = cache;
    do { out.push_back((u8)(c + carry)); c = 0xff; } while (--pending != 0);
    cache = (u8)((u32)low >> 24);
  }
  pending++; low = (u32)low << 8;
}
void RcEnc::bit(u16& p, int b) {
  u32 bound = (range >> RC_BITS) * p;
  if (!b) { range = bound; p += ((1 << RC_BITS) - p) >> RC_MOVE; }
  else { low += bound; range -= bound; p -= p >> RC_MOVE; }
  while (range < RC_TOP) { range <<= 8; shift_low(); }
}
void RcEnc::direct(u32 v, int nb) {
  for (int i = nb - 1; i >= 0; i--) {
    range >>= 1; if ((v >> i) & 1) low += range;
    while (range < RC_TOP) { range <<= 8; shift_low(); }
  }
}
void RcEnc::finish() { for (int i = 0; i < 5; i++) shift_low(); }

void RcDec::init(const u8* s) { p = s + 1; code = 0; range = 0xffffffffu; for (int i = 0; i < 4; i++) code = (code << 8) | *p++; }
int RcDec::bit(u16& pr) {
  u32 bound = (range >> RC_BITS) * pr; int b;
  if (code < bound) { range = bound; pr += ((1 << RC_BITS) - pr) >> RC_MOVE; b = 0; }
  else { code -= bound; range -= bound; pr -= pr >> RC_MOVE; b = 1; }
  while (range < RC_TOP) { range <<= 8; code = (code << 8) | *p++; }
  return b;
}
u32 RcDec::direct(int nb) {
  u32 v = 0;
  for (int i = 0; i < nb; i++) {
    range >>= 1; u32 b = code >= range; if (b) code -= range; v = (v << 1) | b;
    while (range < RC_TOP) { range <<= 8; code = (code << 8) | *p++; }
  }
  return v;
}

// A frequency is stored as its bit length ("class", 0 for an unused symbol) followed by the bits
// below the leading one. The class is coded with a 4-bit tree whose context is the class of the
// previous symbol in this table and a coarse class of the same symbol in the previous table of the
// alphabet; the top remaining bit is modelled per class, the rest are stored plainly.
TabModel::TabModel() {
  for (auto& a : cls) for (auto& b : a) for (auto& c : b) c = 1024;
  for (auto& a : top) for (auto& c : a) c = 1024;
  memset(prev_cls, 0, sizeof prev_cls); have_prev = false;
}
static inline int coarse(int c) { return c == 0 ? 0 : c <= 3 ? 1 : 2; }

// Coarse tables: every frequency except the anchor (the first most frequent symbol) keeps only
// `keep` bits below its leading one; the anchor is not stored and takes whatever is left.
static int table_anchor(const u16* f, int n) { int a = 0; for (int i = 1; i < n; i++) if (f[i] > f[a]) a = i; return a; }
static bool fits_keep(const u16* f, int n, int keep, int a) {
  for (int i = 0; i < n; i++) if (i != a && f[i]) {
    int mb = (int)ilog2(f[i]);
    if (mb > keep && (f[i] & ((1u << (mb - keep)) - 1))) return false;
  }
  return true;
}
static int table_keep(const u16* f, int n, int a) {  // 0 = exact, else 1..3
  for (int k = 1; k <= 3; k++) if (fits_keep(f, n, k, a)) return k;
  return 0;
}

void coarsen_freqs(const u32* cnt, int n, u16* f, int bits) {
  const int SCALE = 1 << bits;
  const int a = table_anchor(f, n);
  if (table_keep(f, n, a)) return;
  double bestgain = 0; int bestk = 0; u16 best[256], g[256];
  for (int k = 1; k <= 3; k++) {
    int sum = 0; double saved = -10, loss = 0;  // 8 bits for the anchor index, 2 for the mode
    for (int i = 0; i < n; i++) {
      g[i] = f[i];
      if (i == a || !f[i]) continue;
      int mb = (int)ilog2(f[i]);
      if (mb > k) { u32 q = 1u << (mb - k); u32 v = (f[i] + q / 2) & ~(q - 1); g[i] = (u16)v; saved += mb - k; }
      sum += g[i];
    }
    int av = SCALE - sum;
    saved += f[a] ? (int)ilog2(f[a]) + 2 : 0;
    if (av < 1) continue;
    g[a] = (u16)av;
    if (table_anchor(g, n) != a || !fits_keep(g, n, k, a)) continue;
    for (int i = 0; i < n; i++) if (cnt[i]) loss += cnt[i] * log2((double)f[i] / g[i]);
    double gain = saved - loss;
    if (gain > bestgain) { bestgain = gain; bestk = k; memcpy(best, g, n * sizeof(u16)); }
  }
  if (bestk) memcpy(f, best, n * sizeof(u16));
}

void write_table(RcEnc& rc, TabModel& tm, const u16* f, int n, int bits) {
  (void)bits;
  const int a = table_anchor(f, n); const int keep = table_keep(f, n, a);
  rc.direct(keep, 2); if (keep) rc.direct(a, 8);
  int pc = 0;
  for (int i = 0; i < n; i++) {
    if (keep && i == a) { tm.prev_cls[i] = 15; continue; }
    int c = f[i] ? (int)ilog2(f[i]) + 1 : 0;
    u16* t = tm.cls[tm.have_prev ? coarse(tm.prev_cls[i]) : 1][pc];
    for (int m = 1, k = 3; k >= 0; k--) { int b = (c >> k) & 1; rc.bit(t[m], b); m = (m << 1) | b; }
    if (c >= 2) {
      int mb = c - 1, kb = keep && mb > keep ? keep : mb;
      u32 r = (f[i] - (1u << mb)) >> (mb - kb);
      rc.bit(tm.top[c][0], (r >> (kb - 1)) & 1);
      if (kb >= 2) rc.direct(r & ((1u << (kb - 1)) - 1), kb - 1);
    }
    tm.prev_cls[i] = (u16)c; pc = c;
  }
  tm.have_prev = true;
}
bool read_table(RcDec& rc, TabModel& tm, u16* f, int n, int bits) {
  const int keep = (int)rc.direct(2); const int a = keep ? (int)rc.direct(8) : -1;
  if (a >= n) return false;
  int pc = 0; u32 sum = 0;
  for (int i = 0; i < n; i++) {
    if (i == a) { tm.prev_cls[i] = 15; f[i] = 0; continue; }
    u16* t = tm.cls[tm.have_prev ? coarse(tm.prev_cls[i]) : 1][pc];
    int m = 1; for (int k = 0; k < 4; k++) m = (m << 1) | rc.bit(t[m]);
    int c = m & 15;
    u32 v = 0;
    if (c == 1) v = 1;
    else if (c >= 2) {
      int mb = c - 1, kb = keep && mb > keep ? keep : mb;
      u32 r = (u32)rc.bit(tm.top[c][0]) << (kb - 1);
      if (kb >= 2) r |= rc.direct(kb - 1);
      v = (1u << mb) + (r << (mb - kb));
    }
    f[i] = (u16)v; sum += v; tm.prev_cls[i] = (u16)c; pc = c;
  }
  if (a >= 0) { if (sum >= (1u << bits)) return false; f[a] = (u16)((1u << bits) - sum); }
  tm.have_prev = true;
  return true;
}

// context map: for each entry a "same as the previous entry" flag, then the id if it differs
void write_map(RcEnc& rc, MapModel& mm, const u8* map, int nctx, int ncl) {
  int mb = 0; while ((1 << mb) < ncl) mb++;
  int prev = 0, ps = 0;
  for (int i = 0; i < nctx; i++) {
    int same = map[i] == prev;
    rc.bit(mm.same[ps], !same);
    if (!same) rc.direct(map[i], mb);
    prev = map[i]; ps = !same;
  }
}
void read_map(RcDec& rc, MapModel& mm, u8* map, int nctx, int ncl) {
  int mb = 0; while ((1 << mb) < ncl) mb++;
  int prev = 0, ps = 0;
  for (int i = 0; i < nctx; i++) {
    int diff = rc.bit(mm.same[ps]);
    if (diff) prev = (int)rc.direct(mb);
    map[i] = (u8)prev; ps = diff;
  }
}

static LitMode g_modes[NUM_LITMODES];
static bool g_modes_init = false;
const LitMode& lit_mode(int m) {
  if (!g_modes_init) {
    memset(g_modes, 0, sizeof g_modes);
    for (int mode = 0; mode < NUM_LITMODES; mode++) {
      LitMode& L = g_modes[mode]; L.nctx = 1; L.sub = 0;
      for (int x = 0; x < 256; x++) {
        switch (mode) {
          case 1: L.A[x] = x; L.nctx = 256; break;                       // previous byte
          case 2: L.B[x] = x; L.nctx = 256; break;                       // match byte
          case 3: L.sub = 1; L.nctx = PMAX; break;                          // position, delta coded
          case 4: L.A[x] = x; L.sub = 1; L.nctx = 256; break;            // previous byte, delta coded
          case 5: L.A[x] = x >> 4; L.B[x] = x & 0xf0; L.nctx = 256; break;
          case 6: L.A[x] = x >> 2; L.nctx = 256; break;                  // position + previous byte
          case 7: L.A[x] = x << 2; L.C[x] = byte_class((u8)x); L.nctx = 1024; break;
          case 8: L.B[x] = x << 2; L.A[x] = byte_class((u8)x); L.nctx = 1024; break;
          case 9: L.B[x] = (x >> 3) << 5; L.A[x] = x >> 3; L.nctx = 1024; break;
          case 10: L.nctx = PMAX; break;                                 // column
          case 11: L.A[x] = x >> 3; L.nctx = PMAX * 32; break;           // column + previous byte
          default: break;
        }
      }
      for (int k = 0; k < PMAX; k++) L.D[k] = mode == 3 || mode == 10 ? k : mode == 11 ? k << 5 : mode == 6 ? (k & 3) << 6 : 0;
    }
    g_modes_init = true;
  }
  return g_modes[m];
}

}  // namespace gt
