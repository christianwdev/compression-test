#include "common.h"
#include <chrono>
using namespace gt;

// Returns the file with PAD zero bytes before and PAD + 64 after the contents.
static std::vector<u8> read_file(const char* path, size_t& n) {
  FILE* f = fopen(path, "rb"); if (!f) { perror(path); exit(1); }
  fseek(f, 0, SEEK_END); n = (size_t)ftell(f); fseek(f, 0, SEEK_SET);
  std::vector<u8> b(n + 2 * PAD + 64, 0);
  if (n && fread(b.data() + PAD, 1, n, f) != n) { perror("read"); exit(1); }
  fclose(f); return b;
}
static void write_file(const char* path, const u8* d, size_t n) {
  FILE* f = fopen(path, "wb"); if (!f) { perror(path); exit(1); }
  if (n && fwrite(d, 1, n, f) != n) { perror("write"); exit(1); }
  fclose(f);
}

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
    else if (!strcmp(a, "-sm")) P.smooth = v; else if (!strcmp(a, "-uc")) P.unseen_cap = atof(argv[i + 1]);
    else if (!strcmp(a, "-f2")) P.far2 = v; else if (!strcmp(a, "-f3")) P.far3 = v; else if (!strcmp(a, "-t")) P.maxtypes = v;
    else { fprintf(stderr, "unknown option %s\n", a); return 2; }
    i++;
  }
  const char* in = argv[argc - 2]; const char* out = argv[argc - 1];
  auto t0 = std::chrono::steady_clock::now();
  if (argv[1][0] == 'c') {
    size_t n = 0; std::vector<u8> srcbuf = read_file(in, n); const u8* src = srcbuf.data() + PAD;
    std::vector<u8> c = compress(src, n, P);
    write_file(out, c.data(), c.size());
    if (P.verbose) fprintf(stderr, "%zu -> %zu (%.3fx) %.2fs\n", n, c.size(), (double)n / c.size(),
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
