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

11. Tried table-based ANS (tANS, three variants) in place of rANS: about 10% slower to decode
    here because the shared bit buffer serialises every symbol. Reverted; kept rANS with five
    states (+4-8% decode).
12. Match-length limit raised to 273 (`nci` -3%, `xml` -1.4%, about 15% more compress time) and
    block headers coded with an adaptive range coder (-13% header bytes).
13. Full corpus (results/gt-v10), clean run against a clean brotli run (results/brotli-11-clean):

    | | brotli -q 11 | gt |
    |---|---|---|
    | Size | 49,564,563 | **48,547,795 (-2.05%)** |
    | Compress | 0.41 MB/s | 0.90 MB/s (2.2x) |
    | Decompress | 233.9 MB/s | 266.2 MB/s (1.14x) |
    | Compress memory (peak) | 238 MB | 221 MB |
    | Decompress memory (peak) | 20 MB | 25 MB |

    Smaller than brotli on 10 of 12 files; `samba` +0.04% and `xml` +0.20% are still larger.
    Decode is faster on 10 of 12; `reymont` (0.95x) and `webster` (0.97x) are slower.

14. Re-tuned the clustering table-cost estimate after headers got cheaper: no gain (lower values
    were 0.07-0.14% worse on the test set). Left unchanged.

15. Coarse frequency tables in headers (fewer significant bits, chosen per table by net saving):
    +0.04% worse on the test set overall (small wins on text, losses elsewhere). Kept in the
    format but off by default.

16. Column period per block: the encoder detects the period P (1..32) at which byte values depend
    most on position mod P and uses position mod P in the contexts of literal-run, match and
    offset symbols and the delta literal mode (previously always pos & 3). The decoder has a
    mask path for P = 4 and a division path otherwise. Two extra column literal modes were
    tried and left disabled (no gain).
    Full corpus (results/gt-v11): **48,233,730 bytes, 2.69% smaller than brotli -q 11**
    (`sao` -7.96%, `mozilla` -1.99%, `samba` -0.98%). Compress 0.89 MB/s (2.2x), decompress
    258 MB/s (1.10x), memory 221 MB / 25 MB. Smaller than brotli on 11 of 12 files; `xml` is
    +0.20%. Decode slower than brotli only on `webster` (0.95x) and `nci` (0.99x).

17. Triple-mapped window ring in the decoder (no slide copy): worse here. Reported memory rose
    from 25 MB to 30 MB and `samba` decoded 15% slower, because each mapping faults separately
    and without huge pages. Reverted.

18. `xml` (0.20% larger than brotli): a fourth parse pass (430,374) or a match-length limit of 512
    (428,576) would beat brotli's 430,566, but each costs 20-40% compress time on every file for
    about 0.01-0.03% overall. Not adopted.

## Ideas not pursued

- x86 call/jump address filter for executables (likely several percent on `mozilla` and
  `ooffice`); left out as file-type specific.
- Learned per-block context quantisers for literals in place of the fixed mode menu.
- Keeping more than one arrival state per position in the optimal parse.

## Open items

- `xml` still 0.20% larger than brotli.
- Decode speed on `webster` (0.95x).
- Decompress memory is 25 MB vs brotli's 20 MB on files larger than the window.
- Compress memory on small files is 10-35% above brotli's (peak is below).
- Literal mode selection still interacts with the parse (up to ~1% on some files).
- Not done: fuzzing the decoder against corrupt input, inputs of 4 GB or more.
