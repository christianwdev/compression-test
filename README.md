# gt

An experimental general-purpose lossless compressor, built to beat brotli at its maximum setting.

## Result

Silesia corpus (12 files, 211,938,580 bytes), same machine, clean runs, every file round-trip
verified by SHA-256 (`results/gt-v11.json` vs `results/brotli-11-clean.json`):

| | brotli 1.1.0 `-q 11` | gt |
|---|---|---|
| Compressed size | 49,564,563 | **48,233,730 (2.69% smaller)** |
| Compress speed | 0.41 MB/s | 0.89 MB/s (2.2x) |
| Decompress speed | 233.9 MB/s | 258.4 MB/s (1.10x) |
| Compress memory (peak) | 238 MB | 221 MB |
| Decompress memory (peak) | 20 MB | 25 MB |

Per file: smaller than brotli on 11 of 12 (`xml` is 0.20% larger); faster to decode on 10 of 12
(`webster` 0.95x, `nci` 0.99x). Compress memory is 10-35% above brotli's on the smaller files.

## Limits

- Research code. The decoder validates table sums and a few header fields but has not been fuzzed;
  do not feed it untrusted input.
- Inputs of 4 GB or more are rejected.
- Tuned and measured on Silesia only.
- Built with clang (faster decoder loop than gcc here); the brotli compared against is the stock
  Ubuntu package.
- The file format changed many times during development and is not stable.

## Build and use

```bash
make
```

```bash
./build/gt c input output.gt
```

```bash
./build/gt d output.gt restored
```

## Reproduce the benchmark

```bash
python3 bench/bench.py --codec brotli-11 --out results/brotli-11.json
```

```bash
python3 bench/bench.py --codec gt --out results/gt.json
```

The corpus is expected at `/home/user/corpus/silesia` (override with `--corpus`).

`GOAL.md` defines the success test. `PROGRESS.md` is the full log of what was tried, including
the ideas that did not work.
