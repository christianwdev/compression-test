# Goal

Build a general-purpose lossless compressor (`gt`) that beats brotli at its maximum setting.

## Success test

Measured by `bench/bench.py` on the Silesia corpus (12 files, 211,938,580 bytes), every result
round-trip verified by SHA-256:

| Requirement | Target (brotli 1.1.0 `-q 11`, this machine) |
|---|---|
| Total compressed size | **<= 49,564,563 bytes** (must) |
| Compression speed | faster than 0.41 MB/s (want) |
| Decompression speed | faster than ~239 MB/s (want) |
| Peak memory | similar: ~238 MB compress, ~20 MB decompress (want) |

"Must" is the ratio. The "want" rows are pursued once ratio is met, without giving ratio back.

## Rules

- No tuning to individual files: one default configuration for all inputs.
- Timings come from the same harness, same machine, same build flags (no `-march=native`).
- Results are only claimed after a verified round-trip.

## How to reproduce

```bash
make
```

```bash
python3 bench/bench.py --codec brotli-11 --out results/brotli-11.json
```

```bash
python3 bench/bench.py --codec gt --out results/gt.json
```

The corpus is expected at `/home/user/corpus/silesia` (override with `--corpus`).

See `PROGRESS.md` for the running log of what has been tried and the current numbers.
