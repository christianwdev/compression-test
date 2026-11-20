#include "common.h"
namespace gt {

void normalize_freqs(const u32* cnt, int n, u16* f) {
  u64 tot = 0; int nz = 0;
  for (int i = 0; i < n; i++) { tot += cnt[i]; nz += cnt[i] != 0; }
  if (tot == 0) { memset(f, 0, n * sizeof(u16)); f[0] = SCALE - 1; f[1] = 1; return; }
  int sum = 0;
  for (int i = 0; i < n; i++) {
    if (!cnt[i]) { f[i] = 0; continue; }
    u64 v = ((u64)cnt[i] * SCALE + tot / 2) / tot;
    if (v < 1) v = 1; if (v > SCALE - 1) v = SCALE - 1;
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

// class (4 bits): 0 = zero, 1..12 = bit length of freq (+ len-1 mantissa bits),
// 13 = zero run 2..9 (3 bits), 14 = zero run 10..265 (8 bits)
void write_table(BitWriter& bw, const u16* f, int n) {
  int i = 0;
  while (i < n) {
    if (f[i] == 0) {
      int run = 1; while (i + run < n && f[i + run] == 0 && run < 265) run++;
      if (run == 1) bw.put(0, 4);
      else if (run < 10) { bw.put(13, 4); bw.put(run - 2, 3); }
      else { bw.put(14, 4); bw.put(run - 10, 8); }
      i += run;
    } else {
      u32 c = ilog2(f[i]) + 1;
      bw.put(c, 4); if (c > 1) bw.put(f[i] - (1u << (c - 1)), c - 1);
      i++;
    }
  }
}
void read_table(BitReader& br, u16* f, int n) {
  int i = 0;
  while (i < n) {
    u32 c = br.get(4);
    if (c == 0) f[i++] = 0;
    else if (c == 13) { int run = br.get(3) + 2; while (run-- && i < n) f[i++] = 0; }
    else if (c == 14) { int run = br.get(8) + 10; while (run-- && i < n) f[i++] = 0; }
    else f[i++] = (u16)((c > 1 ? br.get(c - 1) : 0) + (1u << (c - 1)));
  }
}


// context map: runs of (cluster id, run length); run length = 4-bit bit count + mantissa
void write_map(BitWriter& bw, const u8* map, int nctx, int ncl) {
  int mb = 0; while ((1 << mb) < ncl) mb++;
  int i = 0;
  while (i < nctx) {
    int run = 1; while (i + run < nctx && map[i + run] == map[i]) run++;
    bw.put(map[i], mb);
    u32 c = ilog2(run) + 1; bw.put(c, 4); if (c > 1) bw.put(run - (1u << (c - 1)), c - 1);
    i += run;
  }
}
void read_map(BitReader& br, u8* map, int nctx, int ncl) {
  int mb = 0; while ((1 << mb) < ncl) mb++;
  int i = 0;
  while (i < nctx) {
    u8 id = (u8)br.get(mb); u32 c = br.get(4);
    int run = (int)((c > 1 ? br.get(c - 1) : 0) + (1u << (c - 1)));
    while (run-- && i < nctx) map[i++] = id;
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
          case 3: L.sub = 1; L.nctx = 4; break;                          // position, delta coded
          case 4: L.A[x] = x; L.sub = 1; L.nctx = 256; break;            // previous byte, delta coded
          case 5: L.A[x] = x >> 4; L.B[x] = x & 0xf0; L.nctx = 256; break;
          case 6: L.A[x] = x >> 2; L.nctx = 256; break;                  // position + previous byte
          case 7: L.A[x] = x << 2; L.C[x] = byte_class((u8)x); L.nctx = 1024; break;
          case 8: L.B[x] = x << 2; L.A[x] = byte_class((u8)x); L.nctx = 1024; break;
          case 9: L.B[x] = (x >> 3) << 5; L.A[x] = x >> 3; L.nctx = 1024; break;
          default: break;
        }
      }
      for (int k = 0; k < 4; k++) L.D[k] = mode == 3 ? k : mode == 6 ? k << 6 : 0;
    }
    g_modes_init = true;
  }
  return g_modes[m];
}

}  // namespace gt
