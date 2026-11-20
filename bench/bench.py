#!/usr/bin/env python3
"""Benchmark harness: ratio, compress/decompress time, peak RSS, verified round-trip.

usage: bench.py --codec NAME [--corpus DIR] [--out results/NAME.json] [--files a,b]
Codecs are defined in CODECS as (compress_cmd, decompress_cmd) templates with {in} {out}.
"""
import argparse, json, os, subprocess, sys, tempfile, time, hashlib

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
BIN = os.environ.get("GT_BIN", os.path.join(ROOT, "build", "gt"))

CODECS = {
    "brotli-11":   ("brotli -q 11 -f -o {out} {in}", "brotli -d -f -o {out} {in}"),
    "brotli-11-lw": ("brotli -q 11 --large_window=30 -f -o {out} {in}", "brotli -d -f -o {out} {in}"),
    "brotli-9":    ("brotli -q 9 -f -o {out} {in}", "brotli -d -f -o {out} {in}"),
    "zstd-19":     ("zstd -19 -q -f -o {out} {in}", "zstd -d -q -f -o {out} {in}"),
    "zstd-22":     ("zstd --ultra -22 -q -f -o {out} {in}", "zstd -d -q -f -o {out} {in}"),
    "xz-9":        ("sh -c 'xz -9 -k -c {in} > {out}'", "sh -c 'xz -d -k -c {in} > {out}'"),
}

def run(cmd):
    """Run cmd, return (wall_seconds, peak_rss_kb)."""
    t0 = time.perf_counter()
    p = subprocess.Popen(cmd, shell=True)
    _, status, ru = os.wait4(p.pid, 0)
    t1 = time.perf_counter()
    if status != 0:
        raise RuntimeError(f"command failed ({status}): {cmd}")
    return t1 - t0, ru.ru_maxrss

def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--codec", required=True)
    ap.add_argument("--corpus", default="/home/user/corpus/silesia")
    ap.add_argument("--out")
    ap.add_argument("--files")
    ap.add_argument("--cmd", nargs=2, metavar=("COMP", "DECOMP"),
                    help="ad-hoc codec templates, overrides CODECS lookup")
    ap.add_argument("--dreps", type=int, default=3, help="decompress repetitions (best-of)")
    a = ap.parse_args()
    if a.cmd:
        comp, decomp = a.cmd
    elif a.codec in CODECS:
        comp, decomp = CODECS[a.codec]
    else:
        # gt:<args> -> our codec with extra compress args
        assert a.codec.startswith("gt"), a.codec
        extra = a.codec[3:].replace(",", " ") if len(a.codec) > 3 else ""
        comp = f"{BIN} c {extra} {{in}} {{out}}"
        decomp = f"{BIN} d {{in}} {{out}}"
    files = sorted(os.listdir(a.corpus))
    if a.files:
        want = a.files.split(",")
        files = [f for f in files if f in want]
    rows = []
    tot_in = tot_out = 0
    tot_ct = tot_dt = 0.0
    with tempfile.TemporaryDirectory() as td:
        for f in files:
            src = os.path.join(a.corpus, f)
            cpath = os.path.join(td, f + ".c")
            dpath = os.path.join(td, f + ".d")
            ct, cmem = run(comp.format(**{"in": src, "out": cpath}))
            dt, dmem = 1e9, 0
            for _ in range(a.dreps):
                t, m = run(decomp.format(**{"in": cpath, "out": dpath}))
                dt = min(dt, t); dmem = max(dmem, m)
            ok = sha(src) == sha(dpath)
            n, c = os.path.getsize(src), os.path.getsize(cpath)
            rows.append(dict(file=f, size=n, csize=c, ratio=n / c, ctime=ct, dtime=dt,
                             cmem_mb=cmem / 1024, dmem_mb=dmem / 1024, ok=ok))
            tot_in += n; tot_out += c; tot_ct += ct; tot_dt += dt
            print(f"{f:10s} {n:10d} -> {c:10d}  {n/c:6.3f}x  c {n/ct/1e6:7.2f} MB/s  d {n/dt/1e6:7.1f} MB/s"
                  f"  cmem {cmem/1024:6.0f}M dmem {dmem/1024:5.0f}M  {'OK' if ok else 'ROUNDTRIP FAIL'}", flush=True)
            os.remove(cpath); os.remove(dpath)
            if not ok:
                sys.exit(1)
    summ = dict(codec=a.codec, size=tot_in, csize=tot_out, ratio=tot_in / tot_out,
                c_mbps=tot_in / tot_ct / 1e6, d_mbps=tot_in / tot_dt / 1e6,
                cmem_mb=max(r["cmem_mb"] for r in rows), dmem_mb=max(r["dmem_mb"] for r in rows))
    print(f"TOTAL      {tot_in:10d} -> {tot_out:10d}  {tot_in/tot_out:6.3f}x  c {summ['c_mbps']:7.2f} MB/s  d {summ['d_mbps']:7.1f} MB/s"
          f"  cmem {summ['cmem_mb']:.0f}M dmem {summ['dmem_mb']:.0f}M")
    if a.out:
        os.makedirs(os.path.dirname(a.out), exist_ok=True)
        json.dump(dict(summary=summ, rows=rows), open(a.out, "w"), indent=1)

if __name__ == "__main__":
    main()
