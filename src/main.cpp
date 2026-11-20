#include "common.h"
#include <chrono>
using namespace gt;

int main(int argc, char** argv) {

  if (argc < 4) { fprintf(stderr, "usage: gt c [-w lgwin] [-b lgblock] [-d depth] [-n nice] [-p passes] [-c maxcl] [-v] in out\n       gt d in out\n"); return 2; }
  Params P; int i = 2;
  for (; i < argc - 2; i++) {
    const char* a = argv[i];
    if (!strcmp(a, "-v")) { P.verbose = 1; continue; }
    if (!strcmp(a, "-vv")) { P.verbose = 2; continue; }
    if (i + 1 >= argc - 2) break;
    int v = atoi(argv[i + 1]);
    if (!strcmp(a, "-w")) P.lgwin = v; else if (!strcmp(a, "-b")) P.lgblock = v;
    else if (!strcmp(a, "-d")) P.depth = v; else if (!strcmp(a, "-n")) P.nice = v;
    else if (!strcmp(a, "-p")) P.passes = v; else if (!strcmp(a, "-c")) P.maxcl = v;
    else if (!strcmp(a, "-h")) P.hashbits = v; else if (!strcmp(a, "-m")) P.modemask = v;
    else if (!strcmp(a, "-k")) P.lgchunk = v; else if (!strcmp(a, "-s")) P.lgseg = v; else if (!strcmp(a, "-u")) P.lgsub = v; else if (!strcmp(a, "-rs")) P.rankslots = v;
    else if (!strcmp(a, "-hn")) P.hdr_nz = atof(argv[i + 1]); else if (!strcmp(a, "-hb")) P.hdr_base = atof(argv[i + 1]);
    else if (!strcmp(a, "-ag")) P.align_gain = atof(argv[i + 1]);
    else if (!strcmp(a, "-sm")) P.smooth = v; else if (!strcmp(a, "-uc")) P.unseen_cap = atof(argv[i + 1]);
    else if (!strcmp(a, "-f2")) P.far2 = v; else if (!strcmp(a, "-f3")) P.far3 = v; else if (!strcmp(a, "-t")) P.maxtypes = v;
    else { fprintf(stderr, "unknown option %s\n", a); return 2; }
    i++;
  }
  const char* in = argv[argc - 2]; const char* out = argv[argc - 1];
  if (argv[1][0] == 'b') {  // decode benchmark: gt b <reps> file.gt
    int reps = atoi(out); double best = 1e9;
    for (int r = 0; r < reps; r++) {
      FILE* fi = fopen(in, "rb"); FILE* fo = fopen("/dev/null", "wb"); if (!fi || !fo) return 1;
      setvbuf(fo, nullptr, _IONBF, 0);
      auto t = std::chrono::steady_clock::now();
      if (!decompress(fi, fo)) { fprintf(stderr, "decode error\n"); return 1; }
      double d = std::chrono::duration<double>(std::chrono::steady_clock::now() - t).count();
      if (d < best) best = d;
      fclose(fi); fclose(fo);
    }
    printf("%.4f\n", best); return 0;
  }
  auto t0 = std::chrono::steady_clock::now();
  if (argv[1][0] == 'c') {
    FILE* fi = fopen(in, "rb"); if (!fi) { perror(in); return 1; }
    fseek(fi, 0, SEEK_END); u64 n = (u64)ftell(fi); fseek(fi, 0, SEEK_SET);
    setvbuf(fi, nullptr, _IONBF, 0);
    FILE* fo = fopen(out, "wb"); if (!fo) { perror(out); return 1; }
    u64 csize = 0;
    bool ok = compress(fi, n, fo, P, &csize);
    fclose(fi); if (fclose(fo) != 0) ok = false;
    if (!ok) { fprintf(stderr, "compression failed\n"); return 1; }
    if (P.verbose) fprintf(stderr, "%llu -> %llu (%.3fx) %.2fs\n", (unsigned long long)n, (unsigned long long)csize, (double)n / (csize ? csize : 1),
                           std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count());
  } else {
    FILE* fi = fopen(in, "rb"); if (!fi) { perror(in); return 1; }
    FILE* fo = fopen(out, "wb"); if (!fo) { perror(out); return 1; }
    setvbuf(fo, nullptr, _IONBF, 0);
    bool ok = decompress(fi, fo);
    fclose(fi); if (fclose(fo) != 0) ok = false;
    if (!ok) { fprintf(stderr, "decode error\n"); return 1; }
    if (P.verbose) fprintf(stderr, "decode %.3fs\n", std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count());
  }
  return 0;
}
