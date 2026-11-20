#include "common.h"
#include <chrono>
using namespace gt;

static std::vector<u8> read_file(const char* path) {
  FILE* f = fopen(path, "rb"); if (!f) { perror(path); exit(1); }
  fseek(f, 0, SEEK_END); long n = ftell(f); fseek(f, 0, SEEK_SET);
  std::vector<u8> b(n + 32, 0);
  if (n && fread(b.data(), 1, n, f) != (size_t)n) { perror("read"); exit(1); }
  fclose(f); b.resize(n); return b;
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
    else if (!strcmp(a, "-k")) P.lgchunk = v; else if (!strcmp(a, "-t")) P.maxtypes = v;
    else { fprintf(stderr, "unknown option %s\n", a); return 2; }
    i++;
  }
  const char* in = argv[argc - 2]; const char* out = argv[argc - 1];
  std::vector<u8> src = read_file(in);
  auto t0 = std::chrono::steady_clock::now();
  if (argv[1][0] == 'c') {
    std::vector<u8> c = compress(src.data(), src.size(), P);
    write_file(out, c.data(), c.size());
    if (P.verbose) fprintf(stderr, "%zu -> %zu (%.3fx) %.2fs\n", src.size(), c.size(), (double)src.size() / c.size(),
                           std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count());
  } else {
    std::vector<u8> d; size_t off = 0, len = 0;
    if (!decompress(src.data(), src.size(), d, off, len)) { fprintf(stderr, "decode error\n"); return 1; }
    if (P.verbose) fprintf(stderr, "decode %.3fs\n", std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count());
    write_file(out, d.data() + off, len);
  }
  return 0;
}
