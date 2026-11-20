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

enum { SCALE_BITS = 12, SCALE = 1 << SCALE_BITS, RANS_L = 1 << 16 };
enum { MIN_MATCH = 2, MAX_MATCH = 65537, PAD = 16 };
enum { LL_SYMS = 56, MLC_SYMS = 54, MLK_SYMS = 4 * MLC_SYMS, OFF_SLOTS = 60, NUM_SPECIAL = 13, OFF_SYMS = OFF_SLOTS + NUM_SPECIAL, ALIGN_SYMS = 16 };
enum { NUM_STREAMS = 4, MAXCTX = 1024, MAXCL = 64, NUM_LITMODES = 10, NUM_REPS = 4 };
// Sequence kinds: 0..2 = rep0..rep2, 3 = new offset (slot coded), 4 = rep3, 5..16 = rep0/rep1 +- 1..3.
// Kinds >= 3 are signalled as kind 3 in the MLK symbol; the OFF symbol then selects slot or special.
static const int NEAR_DELTA[6] = {-1, 1, -2, 2, -3, 3};
enum { MLK_NCTX = 12, NT = 16, MAXCH = 4096, LL_MAXCL = 8, MLK_MAXCL = 16, OFF_MAXCL = 16, AL_MAXCL = 4 };
static inline u32 mlk_ctx(u32 ll) {  // context for the match symbol from the literal run length
  return ll < 8 ? ll : ll < 12 ? 8 : ll < 16 ? 9 : ll < 32 ? 10 : 11;
}

static inline u32 ilog2(u32 v) { return 31 - __builtin_clz(v); }
static inline u32 load32(const u8* p) { u32 v; memcpy(&v, p, 4); return v; }
static inline u32 load16(const u8* p) { u16 v; memcpy(&v, p, 2); return v; }
static inline u64 load64(const u8* p) { u64 v; memcpy(&v, p, 8); return v; }

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
static const int LIT_NCTX[NUM_LITMODES] = {1, 256, 256, 4, 256, 256, 256, 1024, 1024, 1024};
static inline bool lit_sub(int mode) { return mode == 3 || mode == 4; }
static inline u32 byte_class(u8 c) {  // coarse class of a byte, used for second-order contexts
  if (c >= 0x80) return 3;
  if ((c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z')) return 2;
  if (c >= 0x30) return 1;
  return 0;
}
static inline u32 lit_ctx(int mode, u32 pos, u8 p1, u8 p2, u8 mb) {
  switch (mode) {
    case 0: return 0;
    case 1: return p1;
    case 2: return mb;
    case 3: return pos & 3;
    case 4: return p1;
    case 5: return (p1 >> 4) | (mb & 0xf0);
    case 6: return ((pos & 3) << 6) | (p1 >> 2);
    case 7: return ((u32)p1 << 2) | byte_class(p2);
    case 8: return ((u32)mb << 2) | byte_class(p1);
    default: return ((u32)(mb >> 3) << 5) | (p1 >> 3);
  }
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
void normalize_freqs(const u32* cnt, int n, u16* f);
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
  int fullevery = 0;   // 1: re-evaluate all literal modes on every model build
  int modetop = 99;     // literal modes fully evaluated after the cheap ranking (>= NUM_LITMODES: all)
  double alpha = 0.1;  // smoothing constant of the cheap ranking
  int lgsub = 20;      // match-storage / price-refinement unit inside a block
  int lgseg = 18;      // optimal-parse segment size
  int far2 = 1 << 8;   // max offset for which a new length-2 match is considered
  int far3 = 1 << 14;  // same for length 3
  int lgchunk = 13;    // table-switch granularity (bytes = 1<<lgchunk)
  int maxtypes = 16;   // max table sets per block for each of literals / commands / distances
  int verbose = 0;
};

// `base` must have PAD readable zero bytes before it and PAD + 8 readable bytes after base + n.
std::vector<u8> compress(const u8* base, size_t n, const Params& p);
// Decodes into buf; the output is buf[out_off .. out_off + out_len).
bool decompress(const u8* in, size_t n, std::vector<u8>& buf, size_t& out_off, size_t& out_len);

}  // namespace gt
