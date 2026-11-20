# Progress log

Baselines on this machine (Silesia total, bytes):

| Codec | Size | Compress MB/s | Decompress MB/s | Compress mem | Decompress mem |
|---|---|---|---|---|---|
| brotli -q 11 | 49,564,563 | 0.41 | 239 | 238 MB | 20 MB |
| xz -9 | 48,795,480 | 1.70 | 59 | 506 MB | 51 MB |
| zstd -19 | 53,024,573 | 1.83 | 518 | 144 MB | 19 MB |

## Design

- LZ77 with a binary-tree match finder (16 MB window) and an iterated optimal parse driven by
  real entropy-coder prices.
- Sequences: literal run, match length joined with offset kind, offset. Four repeat offsets plus
  "repeat offset +-1..3" codes.
- Literals coded under a per-block context mode (previous byte, two-byte, match byte, position,
  subtract-match-byte), with contexts clustered into shared tables.
- Blocks are cut into chunks; each chunk picks a literal / command / distance table set, so
  statistics can change every few kilobytes without paying for new tables.
- Static rANS (12-bit) with 4 interleaved states plus a raw bit stream.

## Log

1. v0: first working codec. Total 50,164,619 (1.2% behind brotli). 1.14 MB/s compress.
2. Instrumented brotli to compare per-component costs. Found its offsets were ~1.5 bits cheaper
   on `xml` thanks to "last distance +-k" codes. Added those plus a 4th repeat offset.
3. Added two-byte literal contexts and clustered match contexts: small gains.
4. Found the main remaining gap was local adaptivity (brotli switches tables inside a block for
   almost no header cost). Added chunk-level table switching: `xml` 460K -> 432K at high settings
   (brotli 430.5K).

## Open items

- Decoder speed (currently ~150 MB/s; not yet optimised).
- Decoder memory: currently holds the whole output; needs a sliding window to match brotli.
- Table header coding is naive (about 2x brotli's header cost).
