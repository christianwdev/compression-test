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

5. Full corpus with `-b 22 -p 3 -t 16` (results/gt-v2a): **48,913,906 bytes, 1.3% smaller than
   brotli -q 11**, 0.86 MB/s compress (2.1x brotli). Not yet acceptable on the other axes:
   657 MB compress memory, 114 MB decompress memory, 140 MB/s decompress.

6. Memory work: compact match storage, segmented parse, 1 MB sub-blocks for match storage,
   streaming encoder (sliding input window) and streaming decoder (window + quarter-window
   headroom). Compress memory 657 MB -> 220 MB; decompress memory 114 MB -> 25 MB.
7. Decoder speed: per-alphabet table precision (11/10/8 bits instead of 12, no ratio loss),
   branch-free rANS refill, raw low offset bits when modelling them does not pay, huge-page hint,
   static link, and a fix for hardened memcpy turning 8-byte stores into library calls (17% of all
   instructions). In-process decode roughly doubled.
8. Literal modes became per chunk type (lookup-table based contexts). Tried and dropped: sampled
   trial parses per mode, smoothed prices for unseen symbols, a tunable table-cost model. Mode
   choice still interacts with the parse and costs up to ~1% on some files (`mr`); open.
9. Position-alignment contexts (pos & 3) for the literal-run, match and offset symbols, chosen per
   block against the position-free alternative. -1.1% on the test set, `sao` -3.9%, text unchanged.
10. Full corpus (results/gt-v8): **48,658,032 bytes, 1.83% smaller than brotli -q 11** and smaller
    than xz -9. 0.89 MB/s compress (2.2x brotli), 220 MB compress memory (brotli 238 MB), 25 MB
    decompress memory (brotli 20 MB). Decode speed about level with brotli; being re-measured.

## Open items

- Literal mode selection is unstable (see 8).
- Longer optimal-parse match limit helps `nci` by ~3% but costs time; needs a cheaper long-match path.
- Table header coding is still naive.
- Decompress memory is 25 MB vs brotli's 20 MB.
