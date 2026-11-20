// gt: experimental general-purpose LZ + context-clustered rANS compressor.
#pragma once
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <vector>
#include <algorithm>

typedef uint8_t u8; typedef uint16_t u16; typedef uint32_t u32; typedef uint64_t u64;

namespace gt {

// Probability precision per alphabet (bits). Smaller tables stay in the fastest cache.
enum { LIT_BITS = 11, MLK_BITS = 11, LL_BITS = 10, OFF_BITS = 10, AL_BITS = 8, MAX_BITS = 12, RANS_L = 1 << 16 };
enum { MIN_MATCH = 2, MAX_MATCH = 65537, PAD = 16 };
enum { LL_SYMS = 56, MLC_SYMS = 54, MLK_SYMS = 4 * MLC_SYMS, OFF_SLOTS = 60, NUM_SPECIAL = 13, OFF_SYMS = OFF_SLOTS + NUM_SPECIAL, ALIGN_SYMS = 16 };
enum { NUM_STREAMS = 4, MAXCTX = 1024, MAXCL = 64, NUM_LITMODES = 10, NUM_REPS = 4 };
// Sequence kinds: 0..2 = rep0..rep2, 3 = new offset (slot coded), 4 = rep3, 5..16 = rep0/rep1 +- 1..3.
// Kinds >= 3 are signalled as kind 3 in the MLK symbol; the OFF symbol then selects slot or special.
static const int NEAR_DELTA[6] = {-1, 1, -2, 2, -3, 3};
enum { MLK_NCTX = 12, NT = 16, MAXCH = 4096, LL_MAXCL = 8, MLK_MAXCL = 16, OFF_MAXCL = 16, AL_MAXCL = 4 };
static const u8 MLK_CTX_LUT[32] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 8, 8, 8, 9, 9, 9, 9, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10};
static inline u32 mlk_ctx(u32 ll) {  // context for the match symbol from the literal run length
  return ll < 32 ? MLK_CTX_LUT[ll] : 11;
}

static inline u32 ilog2(u32 v) { return 31 - __builtin_clz(v); }
static inline u32 load32(const u8* p) { u32 v; __builtin_memcpy(&v, p, 4); return v; }
static inline u32 load16(const u8* p) { u16 v; __builtin_memcpy(&v, p, 2); return v; }
static inline u64 load64(const u8* p) { u64 v; __builtin_memcpy(&v, p, 8); return v; }

// ---- value <-> (symbol, extra bits) codes
static inline void ll_code(u32 v, u32& sym, u32& nb, u32& ex) {
  if (v < 16) { sym = v; nb = 0; ex = 0; return; }
  u32 n = ilog2(v);
  sym = 16 + (n - 4) * 2 + ((v >> (n - 1)) & 1); nb = n - 1; ex = v & ((1u << (n - 1)) - 1);
}
static inline void ml_code(u32 v, u32& sym, u32& nb, u32& ex) {  // v = ml - MIN_MATCH
  if (v < 32) { sym = v; nb = 0; ex = 0; return; }
  u32 n = ilog2(v);
  sym = 32 + (n - 5) * 2 + ((v >> (n - 1)) & 1); nb = n - 1; ex = v & ((1u << (n - 1)) - 1);
}
static inline void off_code(u32 v, u32& sym, u32& nb, u32& ex) {  // v = dist - 1
  if (v < 4) { sym = v; nb = 0; ex = 0; return; }
  u32 n = ilog2(v);
  sym = 2 * n + ((v >> (n - 1)) & 1); nb = n - 1; ex = v & ((1u << (n - 1)) - 1);
}

// ---- literal context modes
// The context of a literal is A[prev byte] + B[match byte] + C[byte before prev] + D[pos & 3], where
// the match byte is the byte at the most recent offset. A mode is a choice of those four tables plus
// whether the literal is coded as (literal - match byte). Each literal type of a block has a mode.
struct LitMode { u16 A[256], B[256], C[256], D[4]; u16 nctx; u8 sub; };
const LitMode& lit_mode(int m);
static inline u32 byte_class(u8 c) {  // coarse class of a byte, used for second-order contexts
  if (c >= 0x80) return 3;
  if ((c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z')) return 2;
  if (c >= 0x30) return 1;
  return 0;
}

// ---- bit IO (LSB first)
struct BitWriter {
  std::vector<u8> buf; u64 acc = 0; int n = 0;
  void put(u32 v, int nb) {
    acc |= (u64)v << n; n += nb;
    while (n >= 8) { buf.push_back((u8)acc); acc >>= 8; n -= 8; }
  }
  void flush() { if (n > 0) { buf.push_back((u8)acc); acc = 0; n = 0; } }
};
struct BitReader {
  const u8* p; u64 acc; int n;
  void init(const u8* s) { p = s; acc = 0; n = 0; }
  inline void refill() { acc |= load64(p) << n; p += (63 - n) >> 3; n |= 56; }
  inline u32 get(int nb) {  // nb <= 32, may be 0
    refill(); u32 v = (u32)(acc & ((1ull << nb) - 1)); acc >>= nb; n -= nb; return v;
  }
};
static inline void put_varint(std::vector<u8>& o, u64 v) {
  while (v >= 128) { o.push_back((u8)(v | 128)); v >>= 7; }
  o.push_back((u8)v);
}
static inline u64 get_varint(const u8*& p) {
  u64 v = 0; int s = 0;
  while (*p & 128) { v |= (u64)(*p++ & 127) << s; s += 7; }
  v |= (u64)(*p++) << s; return v;
}

// ---- frequency tables
void normalize_freqs(const u32* cnt, int n, u16* f, int bits);
void write_table(BitWriter& bw, const u16* f, int n);
void read_table(BitReader& br, u16* f, int n);
void write_map(BitWriter& bw, const u8* map, int nctx, int ncl);
void read_map(BitReader& br, u8* map, int nctx, int ncl);

struct Params {
  int lgwin = 24;      // window = 1<<lgwin
  int lgblock = 21;    // block size
  int depth = 64;      // matchfinder tree depth
  int nice = 128;      // "good enough" match length
  int passes = 3;      // optimal-parse iterations per block
  int maxcl = 32;      // max literal clusters
  int hashbits = 20;
  int modemask = 0xffff; // enabled literal modes
  double hdr_nz = 8.0, hdr_base = 32.0;  // table cost model used when clustering: base + per used symbol
  double align_gain = 0.985;  // model low offset bits only if that costs less than this share of raw
  int smooth = 0;      // smoothed prices for symbols a table has not seen
  double unseen_cap = 14.0;
  int rankslots = 32;  // contexts kept separate when ranking literal modes for a type
  int lgsub = 20;      // match-storage / price-refinement unit inside a block
  int lgseg = 18;      // optimal-parse segment size
  int far2 = 1 << 8;   // max offset for which a new length-2 match is considered
  int far3 = 1 << 14;  // same for length 3
  int lgchunk = 13;    // table-switch granularity (bytes = 1<<lgchunk)
  int maxtypes = 16;   // max table sets per block for each of literals / commands / distances
  int verbose = 0;
};

// Streaming encode of the n bytes of fin to fout. Inputs must be smaller than 4 GB.
bool compress(FILE* fin, u64 n, FILE* fout, const Params& p, u64* out_size);
// Streaming decode from fin to fout.
bool decompress(FILE* fin, FILE* fout);

}  // namespace gt
