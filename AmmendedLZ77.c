/*============================================================================
 *  best_program_nobwt.c  --  C99 port of best_program_nobwt.py, which is
 *      best_program.py ("LZ77 / BWT + context-rANS compressor, v6") with the
 *      BWT pipeline (mode 3) removed: only the stored, entropy-only and
 *      two-domain LZ77 candidates remain.
 *
 *  A transcription, not a redesign: every mechanism, constant, candidate
 *  order and tie-break of the Python is reproduced so that the compressed
 *  output is BYTE-IDENTICAL to best_program_nobwt.py's compress().  Points that
 *  needed care to stay exact:
 *    - Python's sorts are stable; every sort here is stable with the same
 *      keys (quantize, context assignment order, mode shortlist).
 *    - Float sums (_ideal_estimate, stride probe) are accumulated in the
 *      same order as the Python dict iteration (ascending symbol for the
 *      numpy-counted context tables; first appearance for the probe).
 *    - COST[] uses rint(), i.e. round-half-to-even like Python's round().
 *    - LZ hash chains are keyed on the EXACT 4 bytes (Python dict), not on
 *      a lossy hash, so chains and the MAX_CHAIN cap behave identically.
 *    - min() over candidates keeps the FIRST smallest, in the same order.
 *    - The adaptive range coder uses 64-bit integers where Python's are
 *      unbounded; all values stay far below 2^63.
 *  Build with a strict ISO mode (-std=c99 implies -ffp-contract=off in GCC)
 *  so no floating-point contraction changes the model-selection sums.
 *
 *  Container / stream format: see the Python docstring (modes 0, 1, 2).
 *
 *  Build: gcc -O2 -std=c99 -Wall -Wextra -pedantic -o best_program_nobwt best_program_nobwt.c -lm
 *  Usage: best_program_nobwt c <in> <out>      compress
 *         best_program_nobwt d <in> <out>      decompress
 *         best_program_nobwt a <file> [...]    ablation of the delta domain (CSV)
 *         best_program_nobwt f <file> [...]    per-match delta vs whole-file filter (CSV)
 *         best_program_nobwt t <file> [...]    compress + verify, print sizes and
 *                                        every candidate (like LAST_INFO)
 *  Public domain / CC0.
 *==========================================================================*/
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define SCALE_BITS 12
#define MSC        (1 << SCALE_BITS)
#define RANS_L     (1u << 23)
#define WINDOW     65536
#define MAX_LEN    258

/*---- EVOLVE-BLOCK constants (as in the Python) ---------------------------*/
#define MIN_MATCH_RAW      4
#define MIN_MATCH_DELTA    6
#define MAX_CHAIN          120
#define NICE_LEN           160
#define LAZY_MARGIN        0.5
#define INSERT_SKIP        4
static const int DELTA_STRIDES[] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 16 };
#define STRIDE_GATE_MARGIN 0.18
#define SUB_BITS           3
#define SMALL_DATA         4096
static const int SMALL_KS[] = { 1, 2, 3, 4 };
#define ADAPT_MAX          80000
#define ADAPT_O1_MIN       512
#define ADAPT_INC          24
static const int CANDIDATE_KS[] = { 1, 2, 3, 4, 6, 8, 12, 16, 24, 32 };
static const int FRUGAL_KS[] = { 1, 4, 8, 16 };
#define KMEANS_ITERS       15
#define MODES_TO_CLUSTER   2
#define PAIR_BUDGET        180000
#define ABSENT             (24 * 256)
#define NELEM(a) ((int)(sizeof(a) / sizeof((a)[0])))

/*---- ABLATION switches (mode 'a' only; defaults = the Python exactly) -------*/
static int g_min_delta   = MIN_MATCH_DELTA; /* D: delta-domain minimum match       */
static int g_raw_streams = 0;               /* C: 1 = token streams stored raw     */
static int g_explicit    = 0;               /* A: 1 = delta offsets sent explicitly */
static size_t g_last_ndelta = 0;            /* delta-domain matches in last parse  */

typedef int (*CtxFn)(int, int, int);
static int m0 (int a, int b, int c){ (void)a; (void)b; (void)c; return 0; }
static int m1 (int a, int b, int c){ (void)b; (void)c; return a; }
static int m2 (int a, int b, int c){ (void)c; return (a << 4) | (b >> 4); }
static int m3 (int a, int b, int c){ (void)c; return (a & 0xF0) | (b & 0x0F); }
static int m4 (int a, int b, int c){ (void)c; return ((a & 0xC0) >> 3) | (b & 0x3F); }
static int m5 (int a, int b, int c){ (void)c; return (a & 0x0F) | (b & 0xF0); }
static int m6 (int a, int b, int c){ (void)c; return ((a & 0x3F) << 2) | (b >> 6); }
static int m7 (int a, int b, int c){ (void)c; return (a & 0xF8) | ((b & 7) << 4); }
static int m8 (int a, int b, int c){ (void)c; return ((a & 3) << 6) | ((b & 252) >> 2); }
static int m9 (int a, int b, int c){ (void)a; (void)b; return c; }
static int m10(int a, int b, int c){ (void)b; return (c << 4) | (a >> 4); }
static int m11(int a, int b, int c){ (void)b; return (c & 0xF0) | (a >> 4); }
static int m12(int a, int b, int c){ return (a & 0xF8) | ((b & 7) << 4) | ((c & 3) << 1); }
static int m13(int a, int b, int c){ (void)c; return (a & 0x7F) | ((b >> 1) & 0x80); }
static const struct { int n_ctx; CtxFn fn; } MODES[] = {
    {1, m0}, {256, m1}, {4096, m2}, {256, m3}, {256, m4}, {256, m5}, {256, m6},
    {256, m7}, {256, m8}, {256, m9}, {4096, m10}, {256, m11}, {256, m12}, {256, m13} };
#define N_MODES NELEM(MODES)

/*---- utilities ------------------------------------------------------------*/
static void *xmalloc(size_t n){ void *p = malloc(n ? n : 1); if(!p){ fprintf(stderr, "oom\n"); exit(2); } return p; }
static void *xcalloc(size_t n, size_t s){ void *p = calloc(n ? n : 1, s ? s : 1); if(!p){ fprintf(stderr, "oom\n"); exit(2); } return p; }
typedef struct { uint8_t *p; size_t len, cap; } Buf;
static void b_init(Buf *b){ b->p = NULL; b->len = b->cap = 0; }
static void b_free(Buf *b){ free(b->p); b_init(b); }
static void b_reserve(Buf *b, size_t need){
    if(b->cap >= need) return;
    size_t c = b->cap ? b->cap : 64; while(c < need) c *= 2;
    uint8_t *q = realloc(b->p, c); if(!q){ fprintf(stderr, "oom\n"); exit(2); }
    b->p = q; b->cap = c;
}
static void b_push(Buf *b, uint8_t v){ b_reserve(b, b->len + 1); b->p[b->len++] = v; }
static void b_add(Buf *b, const uint8_t *s, size_t n){ if(!n) return; b_reserve(b, b->len + n); memcpy(b->p + b->len, s, n); b->len += n; }
static void b_u32(Buf *b, uint32_t v){ for(int i = 0; i < 4; ++i) b_push(b, (uint8_t)(v >> (8 * i))); }
static uint32_t rd32(const uint8_t *p){ return (uint32_t)p[0] | ((uint32_t)p[1] << 8) | ((uint32_t)p[2] << 16) | ((uint32_t)p[3] << 24); }
static int bitlen(uint64_t v){ int b = 0; while(v){ b++; v >>= 1; } return b; }

/*---- bit I/O (MSB first) --------------------------------------------------*/
typedef struct { Buf buf; uint64_t acc; int n; } BitW;
static void bw_init(BitW *w){ b_init(&w->buf); w->acc = 0; w->n = 0; }
static void bw_write(BitW *w, uint64_t v, int nbits){
    w->acc = (w->acc << nbits) | v; w->n += nbits;
    while(w->n >= 8){ w->n -= 8; b_push(&w->buf, (uint8_t)((w->acc >> w->n) & 0xFF)); }
    w->acc &= ((uint64_t)1 << w->n) - 1;
}
static void bw_value(BitW *w, Buf *out){           /* getvalue(), appended */
    b_add(out, w->buf.p, w->buf.len);
    if(w->n) b_push(out, (uint8_t)((w->acc << (8 - w->n)) & 0xFF));
}
typedef struct { const uint8_t *p; size_t len; uint64_t pos; } BitR;
static uint32_t br_read(BitR *r, int nbits){
    uint32_t v = 0;
    for(int i = 0; i < nbits; ++i){
        uint64_t q = r->pos >> 3;
        int bit = q < r->len ? (r->p[q] >> (7 - (r->pos & 7))) & 1 : 0;
        v = (v << 1) | (uint32_t)bit; r->pos++;
    }
    return v;
}

/*---- stable merge sort of int indices by a key comparator ---------------*/
typedef int (*Cmp)(int, int, const void *);
static void msort(int *a, int n, Cmp cmp, const void *ctx){
    if(n < 2) return;
    int *t = xmalloc((size_t)n * sizeof *t);
    for(int w = 1; w < n; w *= 2){
        for(int lo = 0; lo < n; lo += 2 * w){
            int mid = lo + w < n ? lo + w : n, hi = lo + 2 * w < n ? lo + 2 * w : n;
            int i = lo, j = mid, k = lo;
            while(i < mid && j < hi) t[k++] = cmp(a[j], a[i], ctx) < 0 ? a[j++] : a[i++];
            while(i < mid) t[k++] = a[i++];
            while(j < hi) t[k++] = a[j++];
        }
        memcpy(a, t, (size_t)n * sizeof *a);
    }
    free(t);
}

/*---- quantize --------------------------------------------------------------*/
static int cmp_desc_u64(int x, int y, const void *c){
    const uint64_t *v = c; return v[x] > v[y] ? -1 : v[x] < v[y] ? 1 : 0;
}
static void quantize(const uint64_t counts[256], uint32_t f[256]){
    int items[256], ni = 0; uint64_t total = 0;
    memset(f, 0, 256 * sizeof *f);
    for(int s = 0; s < 256; ++s) if(counts[s]){ items[ni++] = s; total += counts[s]; }
    if(!ni) return;
    long long diff = MSC;
    for(int i = 0; i < ni; ++i){
        uint64_t v = counts[items[i]] * MSC / total; f[items[i]] = v > 1 ? (uint32_t)v : 1; diff -= f[items[i]];
    }
    int order[256]; for(int j = 0; j < ni; ++j) order[j] = items[j];
    msort(order, ni, cmp_desc_u64, counts);
    long long i = 0;
    while(diff){
        int s = order[i % ni];
        if(diff > 0){ f[s]++; diff--; }
        else if(f[s] > 1){ f[s]--; diff++; }
        i++;
    }
}

/*---- frequency-table serialization ---------------------------------------*/
static long freq_rle_bits(const uint32_t *f){
    long bits = 0; int i = 0;
    while(i < 256){
        if(!f[i]){ int j = i; while(j < 256 && !f[j] && j - i < 256) j++; bits += 9; i = j; }
        else { bits += 1 + SCALE_BITS; i++; }
    }
    return bits;
}
static long freq_gamma_bits(const uint32_t *f){
    long bits = 0; int i = 0;
    while(i < 256){
        if(!f[i]){ int j = i; while(j < 256 && !f[j] && j - i < 256) j++; bits += 9; i = j; }
        else { int bl = bitlen(f[i] - 1); bits += 1 + 4 + (bl - 1 > 0 ? bl - 1 : 0); i++; }
    }
    return bits;
}
static long freq_bits(const uint32_t *f){ long a = freq_rle_bits(f), b = freq_gamma_bits(f); return 1 + (a < b ? a : b); }
static void write_freq(BitW *w, const uint32_t *f){
    int gamma = freq_gamma_bits(f) <= freq_rle_bits(f), i = 0;
    bw_write(w, (uint64_t)gamma, 1);
    while(i < 256){
        if(!f[i]){
            int j = i; while(j < 256 && !f[j] && j - i < 256) j++;
            bw_write(w, 0, 1); bw_write(w, (uint64_t)(j - i - 1), 8); i = j;
        } else if(!gamma){
            bw_write(w, 1, 1); bw_write(w, f[i] - 1, SCALE_BITS); i++;
        } else {
            uint32_t v = f[i] - 1; int bl = bitlen(v);
            bw_write(w, 1, 1); bw_write(w, (uint64_t)bl, 4);
            if(bl > 1) bw_write(w, v - (1u << (bl - 1)), bl - 1);
            i++;
        }
    }
}
static void read_freq(BitR *r, uint32_t *f){
    int gamma = (int)br_read(r, 1), n = 0;
    uint32_t tmp[512];
    while(n < 256){
        if(br_read(r, 1) == 0){ int run = (int)br_read(r, 8) + 1; for(int j = 0; j < run; ++j) tmp[n++] = 0; }
        else if(!gamma) tmp[n++] = br_read(r, SCALE_BITS) + 1;
        else {
            int bl = (int)br_read(r, 4); uint32_t v;
            if(bl == 0) v = 0; else if(bl == 1) v = 1; else v = (1u << (bl - 1)) | br_read(r, bl - 1);
            tmp[n++] = v + 1;
        }
    }
    memcpy(f, tmp, 256 * sizeof *f);
}
static long cmap_bits(const int *cmap, int n_ctx, int k){
    int b = bitlen((uint64_t)(k - 1)); if(b < 1) b = 1;
    long raw = (long)n_ctx * b, rle = 0; int i = 0;
    while(i < n_ctx){ int j = i; while(j < n_ctx && cmap[j] == cmap[i] && j - i < 256) j++; rle += b + 8; i = j; }
    return 1 + (raw < rle ? raw : rle);
}
static void write_cmap(BitW *w, const int *cmap, int n_ctx, int k){
    int b = bitlen((uint64_t)(k - 1)); if(b < 1) b = 1;
    long pairs = 0; int i = 0;
    while(i < n_ctx){ int j = i; while(j < n_ctx && cmap[j] == cmap[i] && j - i < 256) j++; pairs++; i = j; }
    if(pairs * (b + 8) < (long)n_ctx * b){
        bw_write(w, 1, 1); i = 0;
        while(i < n_ctx){
            int j = i; while(j < n_ctx && cmap[j] == cmap[i] && j - i < 256) j++;
            bw_write(w, (uint64_t)cmap[i], b); bw_write(w, (uint64_t)(j - i - 1), 8); i = j;
        }
    } else {
        bw_write(w, 0, 1);
        for(i = 0; i < n_ctx; ++i) bw_write(w, (uint64_t)cmap[i], b);
    }
}
static void read_cmap(BitR *r, int *cmap, int n_ctx, int k){
    int b = bitlen((uint64_t)(k - 1)); if(b < 1) b = 1;
    if(br_read(r, 1)){
        int n = 0;
        while(n < n_ctx){
            int v = (int)br_read(r, b), run = (int)br_read(r, 8) + 1;
            for(int j = 0; j < run; ++j){ if(n < n_ctx) cmap[n] = v; n++; }
        }
    } else for(int i = 0; i < n_ctx; ++i) cmap[i] = (int)br_read(r, b);
}

/*---- context counting and model search -----------------------------------*/
static uint32_t *count_ctx(const uint8_t *d, size_t n, int mi){
    int nc = MODES[mi].n_ctx; CtxFn fn = MODES[mi].fn;
    uint32_t *c = xcalloc((size_t)nc * 256, sizeof *c);
    int p1 = 0, p2 = 0, p3 = 0;
    for(size_t i = 0; i < n; ++i){ c[(size_t)fn(p1, p2, p3) * 256 + d[i]]++; p3 = p2; p2 = p1; p1 = d[i]; }
    return c;
}
static double ideal_estimate(const uint32_t *c, int nc){
    double bits = 0.0; int n_active = 0;
    for(int x = 0; x < nc; ++x){
        const uint32_t *row = c + (size_t)x * 256; uint64_t tot = 0;
        for(int s = 0; s < 256; ++s) tot += row[s];
        if(!tot) continue;
        n_active++;
        double lt = log2((double)tot);
        for(int s = 0; s < 256; ++s) if(row[s]) bits += (double)row[s] * (lt - log2((double)row[s]));
    }
    double hg = 550.0 + 35.0 * n_active; if(!(hg < 800.0)) hg = 800.0;
    double pen = 0.18 * nc; if(pen < 1.0) pen = 1.0;
    return bits + hg - pen;
}
static int COST[MSC + 1];
static void init_cost(void){
    for(int f = 1; f <= MSC; ++f){ int v = (int)rint((SCALE_BITS - log2((double)f)) * 256); COST[f] = v > 1 ? v : 1; }
}
typedef struct {                      /* active contexts, ascending ctx id  */
    int nact, nc; int *ctx; int *off; uint8_t *sym; uint32_t *cnt; uint64_t *tot;
} Items;
static void items_build(Items *it, const uint32_t *c, int nc){
    it->nc = nc; it->nact = 0; size_t pairs = 0;
    for(int x = 0; x < nc; ++x){ int any = 0; for(int s = 0; s < 256; ++s) if(c[(size_t)x * 256 + s]){ pairs++; any = 1; } it->nact += any; }
    it->ctx = xmalloc((size_t)it->nact * sizeof(int)); it->off = xmalloc((size_t)(it->nact + 1) * sizeof(int));
    it->sym = xmalloc(pairs); it->cnt = xmalloc(pairs * sizeof(uint32_t)); it->tot = xmalloc((size_t)it->nact * sizeof(uint64_t));
    int a = 0; size_t p = 0;
    for(int x = 0; x < nc; ++x){
        size_t start = p; uint64_t t = 0;
        for(int s = 0; s < 256; ++s){ uint32_t v = c[(size_t)x * 256 + s]; if(v){ it->sym[p] = (uint8_t)s; it->cnt[p] = v; t += v; p++; } }
        if(p > start){ it->ctx[a] = x; it->off[a] = (int)start; it->tot[a] = t; a++; }
    }
    it->off[a] = (int)p;
}
static void items_free(Items *it){ free(it->ctx); free(it->off); free(it->sym); free(it->cnt); free(it->tot); }
static void cluster_tables(const Items *it, const int *cmap, int k, uint32_t (*f)[256], int (*ct)[256]){
    uint64_t (*sums)[256] = xcalloc((size_t)k, sizeof *sums);
    for(int a = 0; a < it->nact; ++a){
        uint64_t *row = sums[cmap[it->ctx[a]]];
        for(int p = it->off[a]; p < it->off[a + 1]; ++p) row[it->sym[p]] += it->cnt[p];
    }
    for(int t = 0; t < k; ++t){
        quantize(sums[t], f[t]);
        for(int s = 0; s < 256; ++s) ct[t][s] = f[t][s] ? COST[f[t][s]] : ABSENT;
    }
    free(sums);
}
static int cmp_tot_desc(int x, int y, const void *c){ const Items *it = c; return it->tot[x] > it->tot[y] ? -1 : it->tot[x] < it->tot[y] ? 1 : 0; }
/* -> payload bits (>> 8 as in Python); cmap (n_ctx) and k_eff out; tables in f */
static uint64_t assign_contexts(const Items *it, int k, int iters, int *cmap, int *k_eff, uint32_t (**fout)[256]){
    int *order = xmalloc((size_t)(it->nact ? it->nact : 1) * sizeof *order);
    for(int a = 0; a < it->nact; ++a) order[a] = a;
    msort(order, it->nact, cmp_tot_desc, it);
    for(int x = 0; x < it->nc; ++x) cmap[x] = 0;
    for(int r = 0; r < it->nact; ++r) cmap[it->ctx[order[r]]] = r % k;
    free(order);
    uint32_t (*f)[256] = xmalloc((size_t)k * sizeof *f);
    int (*ct)[256] = xmalloc((size_t)k * sizeof *ct);
    for(int iter = 0; iter < iters; ++iter){
        cluster_tables(it, cmap, k, f, ct);
        int changed = 0;
        for(int a = 0; a < it->nact; ++a){
            int best_t = -1; uint64_t best = 0;
            for(int t = 0; t < k; ++t){
                uint64_t cost = 0;
                for(int p = it->off[a]; p < it->off[a + 1]; ++p) cost += (uint64_t)it->cnt[p] * (uint64_t)ct[t][it->sym[p]];
                if(best_t < 0 || cost < best){ best_t = t; best = cost; }
            }
            if(best_t != cmap[it->ctx[a]]){ cmap[it->ctx[a]] = best_t; changed = 1; }
        }
        if(!changed) break;
    }
    int remap[256], used[256] = { 0 }, ke = 0;
    for(int a = 0; a < it->nact; ++a) used[cmap[it->ctx[a]]] = 1;
    for(int t = 0; t < 256; ++t) remap[t] = used[t] ? ke++ : 0;
    for(int x = 0; x < it->nc; ++x) cmap[x] = remap[cmap[x]];
    cluster_tables(it, cmap, ke, f, ct);
    uint64_t payload = 0;
    for(int a = 0; a < it->nact; ++a){
        const int *c = ct[cmap[it->ctx[a]]];
        for(int p = it->off[a]; p < it->off[a + 1]; ++p) payload += (uint64_t)it->cnt[p] * (uint64_t)c[it->sym[p]];
    }
    free(ct);
    *k_eff = ke; *fout = f;
    return payload >> 8;
}
static uint64_t choose_k(const Items *it, int frugal, int small, int *best_k, int *best_cmap){
    const int *ks; int nks, iters;
    if(small){ ks = SMALL_KS; nks = NELEM(SMALL_KS); iters = KMEANS_ITERS; }
    else if(frugal){ ks = FRUGAL_KS; nks = NELEM(FRUGAL_KS); iters = 3; }
    else { ks = CANDIDATE_KS; nks = NELEM(CANDIDATE_KS); iters = KMEANS_ITERS; }
    int have = 0, worse = 0; uint64_t best = 0;
    int *cmap = xmalloc((size_t)it->nc * sizeof *cmap);
    for(int i = 0; i < nks; ++i){
        int ke; uint32_t (*f)[256];
        uint64_t payload = assign_contexts(it, ks[i], iters, cmap, &ke, &f);
        uint64_t header = 0;
        for(int t = 0; t < ke; ++t) header += (uint64_t)freq_bits(f[t]);
        if(ke > 1) header += (uint64_t)cmap_bits(cmap, it->nc, ke);
        free(f);
        uint64_t total = header + payload;
        if(!have || total < best){ have = 1; best = total; *best_k = ke; memcpy(best_cmap, cmap, (size_t)it->nc * sizeof *cmap); worse = 0; }
        else if(++worse >= 2) break;
    }
    free(cmap);
    return best;
}
static int cmp_est(int x, int y, const void *c){ const double *e = c; return e[x] < e[y] ? -1 : e[x] > e[y] ? 1 : 0; }
static int select_model(const uint8_t *d, size_t n, int *k_out, int **cmap_out){
    int small = n < SMALL_DATA;
    uint32_t *cnt[N_MODES]; double est[N_MODES]; int order[N_MODES];
    for(int mi = 0; mi < N_MODES; ++mi){ cnt[mi] = count_ctx(d, n, mi); est[mi] = ideal_estimate(cnt[mi], MODES[mi].n_ctx); order[mi] = mi; }
    msort(order, N_MODES, cmp_est, est);
    int have = 0, best_mi = 0; uint64_t best = 0;
    for(int r = 0; r < MODES_TO_CLUSTER && r < N_MODES; ++r){
        int mi = order[r], nc = MODES[mi].n_ctx;
        Items it; items_build(&it, cnt[mi], nc);
        size_t pairs = (size_t)it.off[it.nact];
        int k = 1; int *cmap = xmalloc((size_t)nc * sizeof *cmap);
        uint64_t total = choose_k(&it, pairs > PAIR_BUDGET, small, &k, cmap);
        items_free(&it);
        if(!have || total < best){ have = 1; best = total; best_mi = mi; *k_out = k; free(*cmap_out); *cmap_out = cmap; }
        else free(cmap);
    }
    for(int mi = 0; mi < N_MODES; ++mi) free(cnt[mi]);
    return best_mi;
}

/*---- adaptive headerless range coder -------------------------------------*/
#define RC_TOP  (1ULL << 24)
#define RC_BOT  (1ULL << 16)
#define RC_MASK 0xFFFFFFFFULL
#define ADAPT_LIMIT (1 << 14)
typedef struct { uint32_t counts[256], tree[257], total; } AM;
static void am_rebuild(AM *m){
    memset(m->tree, 0, sizeof m->tree);
    for(int i = 0; i < 256; ++i){ int j = i + 1; m->tree[j] += m->counts[i]; int k = j + (j & -j); if(k <= 256) m->tree[k] += m->tree[j]; }
}
static void am_init(AM *m){ for(int i = 0; i < 256; ++i) m->counts[i] = 1; m->total = 256; am_rebuild(m); }
static uint64_t am_cum(const AM *m, int s){ uint64_t t = 0; int i = s; while(i > 0){ t += m->tree[i]; i -= i & -i; } return t; }
static int am_find(const AM *m, int64_t target){
    int idx = 0, bit = 256;
    while(bit){ int nx = idx + bit; if(nx <= 256 && (int64_t)m->tree[nx] <= target){ idx = nx; target -= m->tree[nx]; } bit >>= 1; }
    return idx;
}
static void am_update(AM *m, int s){
    m->counts[s] += ADAPT_INC; m->total += ADAPT_INC;
    for(int i = s + 1; i <= 256; i += i & -i) m->tree[i] += ADAPT_INC;
    if(m->total > ADAPT_LIMIT){
        uint32_t t = 0;
        for(int i = 0; i < 256; ++i){ uint32_t c = (m->counts[i] + 1) >> 1; m->counts[i] = c; t += c; }
        m->total = t; am_rebuild(m);
    }
}
static AM *am_get(AM **tab, int key){ if(!tab[key]){ tab[key] = xmalloc(sizeof(AM)); am_init(tab[key]); } return tab[key]; }
static void adaptive_compress(const uint8_t *d, size_t n, int order1, Buf *out){
    uint64_t low = 0, rng = RC_MASK; AM *tab[256] = { 0 }; int prev = 0;
    for(size_t i = 0; i < n; ++i){
        int b = d[i]; AM *m = am_get(tab, order1 ? prev : 0);
        uint64_t r = rng / m->total;
        low += r * am_cum(m, b); rng = r * m->counts[b];
        for(;;){
            if((low ^ (low + rng)) < RC_TOP) {}
            else if(rng < RC_BOT) rng = (0 - low) & (RC_BOT - 1);
            else break;
            b_push(out, (uint8_t)((low >> 24) & 0xFF));
            low = (low << 8) & RC_MASK; rng = (rng << 8) & RC_MASK;
        }
        am_update(m, b);
        if(order1) prev = b;
    }
    for(int i = 0; i < 4; ++i){ b_push(out, (uint8_t)((low >> 24) & 0xFF)); low = (low << 8) & RC_MASK; }
    for(int i = 0; i < 256; ++i) free(tab[i]);
}
static int64_t floordiv(int64_t a, int64_t b){ int64_t q = a / b; if((a % b != 0) && ((a < 0) != (b < 0))) q--; return q; }
static size_t adaptive_decompress(const uint8_t *blob, size_t end, size_t pos, size_t n, int order1, uint8_t *out){
    int64_t code = 0; size_t p = pos;
    for(int i = 0; i < 4; ++i){ code = (code << 8) | (p < end ? blob[p] : 0); p++; }
    int64_t low = 0, rng = (int64_t)RC_MASK; AM *tab[256] = { 0 }; int prev = 0;
    for(size_t i = 0; i < n; ++i){
        AM *m = am_get(tab, order1 ? prev : 0);
        int64_t tot = m->total, r = rng / tot;
        int64_t target = floordiv(code - low, r);
        if(target >= tot) target = tot - 1;
        int s = am_find(m, target);
        low += r * (int64_t)am_cum(m, s); rng = r * (int64_t)m->counts[s];
        for(;;){
            if((low ^ (low + rng)) < (int64_t)RC_TOP) {}
            else if(rng < (int64_t)RC_BOT) rng = (-low) & (int64_t)(RC_BOT - 1);
            else break;
            code = ((code << 8) | (p < end ? blob[p] : 0)) & (int64_t)RC_MASK; p++;
            low = (low << 8) & (int64_t)RC_MASK; rng = (rng << 8) & (int64_t)RC_MASK;
        }
        am_update(m, s);
        out[i] = (uint8_t)s;
        if(order1) prev = s;
    }
    for(int i = 0; i < 256; ++i) free(tab[i]);
    return p;
}

/*---- entropy engine --------------------------------------------------------*/
static void entropy_compress(const uint8_t *d, size_t n, Buf *out){
    if(n == 0){ b_push(out, 0); return; }
    if(g_raw_streams){ b_push(out, 0); b_add(out, d, n); return; }   /* ablation C */
    int k = 1; int *cmap = NULL;
    int mi = select_model(d, n, &k, &cmap);
    int nc = MODES[mi].n_ctx; CtxFn fn = MODES[mi].fn;
    uint64_t (*sums)[256] = xcalloc((size_t)k, sizeof *sums);
    uint8_t *tid = xmalloc(n);
    { int p1 = 0, p2 = 0, p3 = 0;
      for(size_t i = 0; i < n; ++i){ int t = cmap[fn(p1, p2, p3)]; tid[i] = (uint8_t)t; sums[t][d[i]]++; p3 = p2; p2 = p1; p1 = d[i]; } }
    uint32_t (*fq)[256] = xmalloc((size_t)k * sizeof *fq), (*cum)[256] = xmalloc((size_t)k * sizeof *cum);
    for(int t = 0; t < k; ++t){ quantize(sums[t], fq[t]); uint32_t a = 0; for(int s = 0; s < 256; ++s){ cum[t][s] = a; a += fq[t][s]; } }
    free(sums);
    Buf rv; b_init(&rv);
    uint64_t x = RANS_L;
    for(size_t i = n; i-- > 0; ){
        int t = tid[i], s = d[i]; uint64_t fv = fq[t][s], xmax = fv << 19;
        while(x >= xmax){ b_push(&rv, (uint8_t)(x & 0xFF)); x >>= 8; }
        x = ((x / fv) << SCALE_BITS) + x % fv + cum[t][s];
    }
    BitW w; bw_init(&w);
    bw_write(&w, (uint64_t)mi, 8); bw_write(&w, (uint64_t)k, 8);
    if(k > 1) write_cmap(&w, cmap, nc, k);
    for(int t = 0; t < k; ++t) write_freq(&w, fq[t]);
    Buf full; b_init(&full); b_push(&full, 1); bw_value(&w, &full); b_free(&w.buf);
    b_u32(&full, (uint32_t)x);
    for(size_t i = rv.len; i-- > 0; ) b_push(&full, rv.p[i]);
    b_free(&rv); free(fq); free(cum); free(tid); free(cmap);
    /* candidates in Python order: stored, full, adaptive o0, adaptive o1 */
    Buf best; b_init(&best); b_push(&best, 0); b_add(&best, d, n);
    if(full.len < best.len){ b_free(&best); best = full; } else b_free(&full);
    for(int o = 0; o < 2; ++o){
        if(n > ADAPT_MAX || (o == 1 && n < ADAPT_O1_MIN)) continue;
        Buf a; b_init(&a); b_push(&a, 2); b_push(&a, (uint8_t)o); adaptive_compress(d, n, o, &a);
        if(a.len < best.len){ b_free(&best); best = a; } else b_free(&a);
    }
    b_add(out, best.p, best.len); b_free(&best);
}
static size_t entropy_decompress(const uint8_t *blob, size_t end, size_t pos, size_t n, uint8_t *out){
    if(n == 0) return pos + 1;
    int inner = blob[pos];
    if(inner == 0){ memcpy(out, blob + pos + 1, n); return pos + 1 + n; }
    if(inner == 2) return adaptive_decompress(blob, end, pos + 2, n, blob[pos + 1], out);
    BitR r = { blob, end, (uint64_t)(pos + 1) * 8 };
    int mi = (int)br_read(&r, 8), k = (int)br_read(&r, 8);
    int nc = MODES[mi].n_ctx; CtxFn fn = MODES[mi].fn;
    int *cmap = xcalloc((size_t)nc, sizeof *cmap);
    if(k > 1) read_cmap(&r, cmap, nc, k);
    uint32_t (*fq)[256] = xmalloc((size_t)k * sizeof *fq), (*cum)[256] = xmalloc((size_t)k * sizeof *cum);
    uint8_t *s2s = xmalloc((size_t)k * MSC);
    for(int t = 0; t < k; ++t){
        read_freq(&r, fq[t]); uint32_t a = 0;
        for(int s = 0; s < 256; ++s){ cum[t][s] = a; for(uint32_t j = 0; j < fq[t][s]; ++j) s2s[(size_t)t * MSC + a + j] = (uint8_t)s; a += fq[t][s]; }
    }
    size_t base = (size_t)((r.pos + 7) >> 3);
    uint64_t x = rd32(blob + base); size_t p = base + 4;
    int p1 = 0, p2 = 0, p3 = 0;
    for(size_t i = 0; i < n; ++i){
        int t = cmap[fn(p1, p2, p3)];
        uint32_t slot = (uint32_t)(x & (MSC - 1));
        int s = s2s[(size_t)t * MSC + slot];
        x = (uint64_t)fq[t][s] * (x >> SCALE_BITS) + slot - cum[t][s];
        while(x < RANS_L){ x = (x << 8) | (p < end ? blob[p] : 0); p++; }
        out[i] = (uint8_t)s; p3 = p2; p2 = p1; p1 = s;
    }
    free(cmap); free(fq); free(cum); free(s2s);
    return p;
}

/*---- stride probe ------------------------------------------------------------*/
static double sample_entropy(const uint8_t *d, size_t n, int stride, size_t step){
    uint64_t cnt[256] = { 0 }; int ord[256], no = 0; uint64_t total = 0;
    for(size_t i = (size_t)stride; i < n; i += step){
        int v = stride ? (d[i] - d[i - (size_t)stride]) & 0xFF : d[i];
        if(!cnt[v]) ord[no++] = v;
        cnt[v]++; total++;
    }
    double h = 0.0, lt = log2((double)total);
    for(int j = 0; j < no; ++j) h += (double)cnt[ord[j]] * (lt - log2((double)cnt[ord[j]]));
    return h / (double)total;
}
static int choose_stride(const uint8_t *d, size_t n){
    size_t lim = n < 256 ? n : 256, step = n / 65536; if(step < 1) step = 1;
    int any = 0; for(int i = 0; i < NELEM(DELTA_STRIDES); ++i) if(DELTA_STRIDES[i] > 0 && (size_t)DELTA_STRIDES[i] < lim) any = 1;
    if(!any) return 0;
    double raw_h = sample_entropy(d, n, 0, step), best_h = 0; int best_s = 0, have = 0;
    for(int i = 0; i < NELEM(DELTA_STRIDES); ++i){
        int s = DELTA_STRIDES[i];
        if(!(s > 0 && (size_t)s < lim)) continue;
        double h = sample_entropy(d, n, s, step);
        if(!have || h < best_h){ have = 1; best_s = s; best_h = h; }
    }
    return (have && best_h < raw_h - STRIDE_GATE_MARGIN) ? best_s : 0;
}

/*---- exact-key hash chains (the Python dict of 4-byte keys) --------------*/
typedef struct { uint32_t *key; int32_t *pos; size_t cap, used; } Map;
static uint32_t mix32(uint32_t k){ k ^= k >> 16; k *= 0x7feb352dU; k ^= k >> 15; k *= 0x846ca68bU; k ^= k >> 16; return k; }
static void map_init(Map *m, size_t cap){ m->cap = cap; m->used = 0; m->key = xmalloc(cap * sizeof *m->key); m->pos = xmalloc(cap * sizeof *m->pos); for(size_t i = 0; i < cap; ++i) m->pos[i] = -1; }
static void map_free(Map *m){ free(m->key); free(m->pos); }
static size_t map_slot(const Map *m, uint32_t k){ size_t i = mix32(k) & (m->cap - 1); while(m->pos[i] >= 0 && m->key[i] != k) i = (i + 1) & (m->cap - 1); return i; }
static void map_put(Map *m, uint32_t k, int32_t p);
static void map_grow(Map *m){ Map g; map_init(&g, m->cap * 2); for(size_t i = 0; i < m->cap; ++i) if(m->pos[i] >= 0) map_put(&g, m->key[i], m->pos[i]); map_free(m); *m = g; }
static void map_put(Map *m, uint32_t k, int32_t p){
    size_t i = map_slot(m, k);
    if(m->pos[i] < 0){ if(2 * (m->used + 1) > m->cap){ map_grow(m); i = map_slot(m, k); } m->used++; m->key[i] = k; }
    m->pos[i] = p;
}
typedef struct { const uint8_t *buf; size_t n; int minm; Map head; int32_t *prev; } Finder;
static uint32_t k4(const uint8_t *p){ return (uint32_t)p[0] | ((uint32_t)p[1] << 8) | ((uint32_t)p[2] << 16) | ((uint32_t)p[3] << 24); }
static void fd_init(Finder *f, const uint8_t *buf, size_t n, int minm){
    f->buf = buf; f->n = n; f->minm = minm < 4 ? 4 : minm; map_init(&f->head, 1024);
    f->prev = xmalloc((n ? n : 1) * sizeof *f->prev); for(size_t i = 0; i < n; ++i) f->prev[i] = -1;
}
static void fd_free(Finder *f){ map_free(&f->head); free(f->prev); }
static void fd_insert(Finder *f, size_t pos){
    if(pos + 4 <= f->n){ uint32_t k = k4(f->buf + pos); size_t sl = map_slot(&f->head, k); f->prev[pos] = f->head.pos[sl]; map_put(&f->head, k, (int32_t)pos); }
}
static void fd_find(Finder *f, size_t pos, int64_t lo_limit, int *len, int *dist){
    *len = 0; *dist = 0;
    if(pos + 4 > f->n) return;
    const uint8_t *b = f->buf;
    int64_t cand = f->head.pos[map_slot(&f->head, k4(b + pos))];
    int best = 0, bd = 0, chain = 0;
    int64_t limit = (int64_t)pos - WINDOW; if(lo_limit > limit) limit = lo_limit;
    int maxl = pos + MAX_LEN <= f->n ? MAX_LEN : (int)(f->n - pos);
    while(cand >= limit && chain < MAX_CHAIN){
        if(best < maxl && b[cand + best] == b[pos + best]){
            int l = 0; while(l < maxl && b[cand + l] == b[pos + l]) l++;
            if(l > best){ best = l; bd = (int)((int64_t)pos - cand); if(l >= NICE_LEN) break; }
        }
        cand = f->prev[cand]; chain++;
    }
    if(best >= f->minm){ *len = best; *dist = bd; }
}

/*---- two-domain LZ77 -----------------------------------------------------------*/
static void dist_band(int d, int sb, int *sym, uint32_t *ex, int *exn){
    int bl = bitlen((uint64_t)d);
    if(bl - 1 <= sb){ *sym = ((bl - 1) << sb) | (d - (1 << (bl - 1))); *ex = 0; *exn = 0; return; }
    int en = bl - 1 - sb;
    *sym = ((bl - 1) << sb) | ((d >> en) & ((1 << sb) - 1)); *ex = (uint32_t)d & ((1u << en) - 1); *exn = en;
}
static int dist_unband(int sym, BitR *r, int sb){
    int blm1 = sym >> sb, sub = sym & ((1 << sb) - 1);
    if(blm1 <= sb) return (1 << blm1) + sub;
    int en = blm1 - sb;
    return (1 << blm1) | (sub << en) | (int)br_read(r, en);
}
typedef struct { Buf flags, lits, lens, dbands, offs; BitW extras; } Toks;
typedef struct { Finder r, d; int has_d, stride; } Two;
static void best_at(Two *T, size_t pos, int *l, int *d, int *kind){
    fd_find(&T->r, pos, 0, l, d); *kind = 1;
    if(T->has_d && pos > (size_t)T->stride){
        int l2, d2; fd_find(&T->d, pos, T->stride, &l2, &d2);
        if(l2 > *l && (int64_t)pos - d2 - T->stride >= 0){ *l = l2; *d = d2; *kind = 2; }
    }
}
static void ins2(Two *T, size_t pos){ fd_insert(&T->r, pos); if(T->has_d) fd_insert(&T->d, pos); }
static void tokenize(const uint8_t *data, size_t n, int stride, Toks *K){
    uint8_t *diff = NULL; Two T; T.stride = stride; T.has_d = stride != 0;
    b_init(&K->flags); b_init(&K->lits); b_init(&K->lens); b_init(&K->dbands); b_init(&K->offs); bw_init(&K->extras);
    fd_init(&T.r, data, n, MIN_MATCH_RAW);
    if(stride){
        diff = xcalloc(n, 1);
        for(size_t i = (size_t)stride; i < n; ++i) diff[i] = (uint8_t)(data[i] - data[i - (size_t)stride]);
        fd_init(&T.d, diff, n, g_min_delta);
    }
    size_t i = 0, last_match = 0;
    while(i < n){
        int l, d, kind; best_at(&T, i, &l, &d, &kind);
        if(l){
            ins2(&T, i);
            if(i + 1 < n){
                int l2, d2, k2; best_at(&T, i + 1, &l2, &d2, &k2);
                if(l2 > l + LAZY_MARGIN){ b_push(&K->flags, 0); b_push(&K->lits, data[i]); i++; continue; }
            }
            b_push(&K->flags, (uint8_t)kind); b_push(&K->lens, (uint8_t)(l - 4));
            int sym, exn; uint32_t ex; dist_band(d, SUB_BITS, &sym, &ex, &exn);
            b_push(&K->dbands, (uint8_t)sym);
            if(kind == 2)       /* the s offsets a delta-domain match implies (ablation A) */
                for(int j = 0; j < stride; ++j)
                    b_push(&K->offs, (uint8_t)(data[i - (size_t)stride + (size_t)j] - data[i - (size_t)d - (size_t)stride + (size_t)j]));
            if(exn) bw_write(&K->extras, ex, exn);
            size_t step = l > 32 ? INSERT_SKIP : 1;
            for(size_t j = i + 1; j < i + (size_t)l; j += step) ins2(&T, j);
            i += (size_t)l; last_match = i;
        } else {
            b_push(&K->flags, 0); b_push(&K->lits, data[i]);
            if(INSERT_SKIP <= 1 || (i - last_match) < 256 || (i & (INSERT_SKIP - 1)) == 0) ins2(&T, i);
            i++;
        }
    }
    fd_free(&T.r); if(stride){ fd_free(&T.d); free(diff); }
}
static int detokenize(size_t n, int stride, int sb, const uint8_t *flags, size_t nt, const uint8_t *lits,
                      const uint8_t *lens, const uint8_t *dbands, const uint8_t *offs, BitR *er, uint8_t *out){
    size_t li = 0, mi = 0, p = 0, dm = 0;
    for(size_t k = 0; k < nt; ++k){
        if(flags[k] == 0){ if(p >= n) return 0; out[p++] = lits[li++]; }
        else {
            size_t l = (size_t)lens[mi] + 4; size_t d = (size_t)dist_unband(dbands[mi], er, sb); mi++;
            if(d > p || p + l > n) return 0;
            size_t src = p - d;
            if(flags[k] == 1) for(size_t t = 0; t < l; ++t) out[p + t] = out[src + t];
            else {
                size_t s = (size_t)stride; if(src < s || p < s || s == 0) return 0;
                if(offs){                          /* explicit offsets (ablation A) */
                    const uint8_t *o = offs + dm * s;
                    for(size_t t = 0; t < l; ++t) out[p + t] = (uint8_t)(out[src + t] + o[t % s]);
                } else
                    for(size_t t = 0; t < l; ++t) out[p + t] = (uint8_t)(out[p + t - s] + out[src + t] - out[src + t - s]);
                dm++;
            }
            p += l;
        }
    }
    return p == n;
}

/*---- drivers -------------------------------------------------------------------*/
typedef struct { size_t stored, entropy, lz, lz0; int stride, has_lz0, has_lz; const char *chosen; } Info;
static Info g_info;
static void lz_candidate(const uint8_t *d, size_t n, int stride, Buf *out){
    Toks K; tokenize(d, n, stride, &K);
    int expl = g_explicit && stride;       /* sub_bits bit 7 flags the offsets stream */
    b_u32(out, (uint32_t)n); b_push(out, 1); b_push(out, (uint8_t)stride); b_push(out, (uint8_t)(SUB_BITS | (expl ? 0x80 : 0)));
    b_u32(out, (uint32_t)K.flags.len);
    entropy_compress(K.flags.p, K.flags.len, out);
    entropy_compress(K.lits.p, K.lits.len, out);
    entropy_compress(K.lens.p, K.lens.len, out);
    entropy_compress(K.dbands.p, K.dbands.len, out);
    if(expl) entropy_compress(K.offs.p, K.offs.len, out);
    bw_value(&K.extras, out);
    g_last_ndelta = stride ? K.offs.len / (size_t)stride : 0;
    b_free(&K.offs); b_free(&K.flags); b_free(&K.lits); b_free(&K.lens); b_free(&K.dbands); b_free(&K.extras.buf);
}
static void compress(const uint8_t *d, size_t n, Buf *out){
    memset(&g_info, 0, sizeof g_info);
    if(n == 0){ b_u32(out, 0); g_info.chosen = "empty"; return; }
    Buf c[4]; const char *name[4]; int nc = 0;
    b_init(&c[nc]); b_u32(&c[nc], (uint32_t)n); b_push(&c[nc], 0); b_add(&c[nc], d, n); name[nc++] = "stored";
    b_init(&c[nc]); b_u32(&c[nc], (uint32_t)n); b_push(&c[nc], 2); entropy_compress(d, n, &c[nc]); name[nc++] = "entropy";
    g_info.stored = c[0].len; g_info.entropy = c[1].len;
    if(n >= 16){
        int stride = choose_stride(d, n); g_info.stride = stride; g_info.has_lz = 1;
        b_init(&c[nc]); lz_candidate(d, n, stride, &c[nc]); g_info.lz = c[nc].len; name[nc++] = stride ? "lz(stride>0)" : "lz(stride=0)";
        if(stride){ b_init(&c[nc]); lz_candidate(d, n, 0, &c[nc]); g_info.lz0 = c[nc].len; g_info.has_lz0 = 1; name[nc++] = "lz(stride=0)"; }
    }
    int best = 0;                                /* first smallest, as min() */
    for(int i = 1; i < nc; ++i) if(c[i].len < c[best].len) best = i;
    g_info.chosen = name[best];
    b_add(out, c[best].p, c[best].len);
    for(int i = 0; i < nc; ++i) b_free(&c[i]);
}
static int decompress(const uint8_t *blob, size_t len, Buf *out){
    if(len < 4) return 0;
    size_t n = rd32(blob); b_reserve(out, n ? n : 1); out->len = n;
    if(n == 0) return 1;
    if(len < 5) return 0;
    int mode = blob[4];
    if(mode == 0){ if(len < 5 + n) return 0; memcpy(out->p, blob + 5, n); return 1; }
    if(mode == 2){ entropy_decompress(blob, len, 5, n, out->p); return 1; }
    if(mode != 1 || len < 11) return 0;
    size_t pos = 5; int stride = blob[pos++], sb = blob[pos++], expl = sb & 0x80; sb &= 0x7F;
    size_t nt = rd32(blob + pos); pos += 4;
    uint8_t *fl = xmalloc(nt ? nt : 1);
    pos = entropy_decompress(blob, len, pos, nt, fl);
    size_t nm = 0; for(size_t i = 0; i < nt; ++i) nm += fl[i] != 0;
    uint8_t *li = xmalloc(nt - nm + 1), *le = xmalloc(nm + 1), *db = xmalloc(nm + 1);
    pos = entropy_decompress(blob, len, pos, nt - nm, li);
    pos = entropy_decompress(blob, len, pos, nm, le);
    pos = entropy_decompress(blob, len, pos, nm, db);
    uint8_t *of = NULL;
    if(expl){
        size_t nd = 0; for(size_t i = 0; i < nt; ++i) nd += fl[i] == 2;
        of = xmalloc(nd * (size_t)stride + 1);
        pos = entropy_decompress(blob, len, pos, nd * (size_t)stride, of);
    }
    BitR er = { blob, len, (uint64_t)pos * 8 };
    int ok = detokenize(n, stride, sb, fl, nt, li, le, db, of, &er, out->p);
    free(of); free(fl); free(li); free(le); free(db);
    return ok;
}

/*---- CLI -----------------------------------------------------------------------*/
static uint8_t *read_all(const char *path, size_t *n){
    FILE *f = fopen(path, "rb"); if(!f){ perror(path); return NULL; }
    size_t cap = 1 << 16, len = 0, r; uint8_t *b = xmalloc(cap);
    while((r = fread(b + len, 1, cap - len, f)) > 0){ len += r; if(len == cap){ cap *= 2; uint8_t *q = realloc(b, cap); if(!q){ fprintf(stderr, "oom\n"); exit(2); } b = q; } }
    fclose(f); *n = len; return b;
}
static int write_all(const char *path, const uint8_t *p, size_t n){
    FILE *f = fopen(path, "wb"); if(!f){ perror(path); return 0; }
    int ok = n == 0 || fwrite(p, 1, n, f) == n; return fclose(f) == 0 && ok;
}
/*---- ablation of the delta-domain LZ77 (mode 'a') ------------------------------
 * For each file, LZ77-pipeline sizes (container included) with the factors of
 * the delta domain switched off one at a time; every stream is decoded and
 * compared with the input.  Variants use the probed stride P unless noted,
 * and fall back to the stride-0 candidate when the probe disables the delta
 * domain (P = 0), except stride1_bytes, which forces stride 1 on every file.
 *   lz0_bytes       no delta domain (stride 0)
 *   full_bytes      as shipped: stride P, implicit offsets, entropy, min 6
 *   explicit_bytes  A: same parse, the P offsets per delta match are sent
 *   stride1_bytes   B: stride 1 instead of the probed stride, on every file
 *   minD4_bytes     D: delta-domain minimum match 4 instead of 6
 *   raw_lz0_bytes / raw_full_bytes   C: token streams stored, no entropy
 *-----------------------------------------------------------------------------*/
static size_t lz_verified(const uint8_t *d, size_t n, int stride, int *fail){
    Buf c; b_init(&c); lz_candidate(d, n, stride, &c);
    Buf r; b_init(&r);
    if(!(decompress(c.p, c.len, &r) && r.len == n && (n == 0 || memcmp(r.p, d, n) == 0))) *fail = 1;
    size_t len = c.len; b_free(&c); b_free(&r);
    return len;
}
static int ablation(int nf, char **files){
    printf("file,n,probe_stride,lz0_bytes,full_bytes,explicit_bytes,stride1_bytes,minD4_bytes,"
           "raw_lz0_bytes,raw_full_bytes,delta_matches_full,delta_matches_stride1,verified\n");
    int rc = 0;
    for(int i = 0; i < nf; ++i){
        size_t n; uint8_t *d = read_all(files[i], &n);
        if(!d){ rc = 2; continue; }
        if(n < 16){ free(d); continue; }                 /* no LZ candidate below 16 */
        int fail = 0, P = choose_stride(d, n);
        size_t lz0 = lz_verified(d, n, 0, &fail);
        size_t full = P ? lz_verified(d, n, P, &fail) : lz0, ndf = P ? g_last_ndelta : 0;
        g_explicit = 1;  size_t ex = P ? lz_verified(d, n, P, &fail) : lz0;  g_explicit = 0;
        size_t s1 = lz_verified(d, n, 1, &fail), nd1 = g_last_ndelta;
        g_min_delta = 4; size_t m4 = P ? lz_verified(d, n, P, &fail) : lz0; g_min_delta = MIN_MATCH_DELTA;
        g_raw_streams = 1;
        size_t r0 = lz_verified(d, n, 0, &fail), rf = P ? lz_verified(d, n, P, &fail) : r0;
        g_raw_streams = 0;
        const char *b = files[i] + strlen(files[i]);
        while(b > files[i] && b[-1] != '/' && b[-1] != '\\') b--;
        printf("%s,%zu,%d,%zu,%zu,%zu,%zu,%zu,%zu,%zu,%zu,%zu,%s\n", b, n, P, lz0, full, ex, s1, m4, r0, rf, ndf, nd1, fail ? "FAIL" : "ok");
        fflush(stdout);
        if(fail) rc = 4;
        free(d);
    }
    return rc;
}

/*---- per-match difference domain vs whole-file delta filter (mode 'f') -------
 * Same LZ77 + entropy back end for both.  The whole-file filter is the classic
 * delta filter (xz --delta, PNG "Sub"): y[i] = x[i] - x[i-s] for i >= s, first
 * s bytes kept; it is then coded by the stride-0 LZ77 pipeline (and, for the
 * full-system columns, also by the entropy-only mode) and charged +1 byte for
 * signalling s.  "auto" uses the stride the probe picks (P); "oracle" takes
 * the best of all DELTA_STRIDES in hindsight.  Every column keeps the no-delta
 * candidate as a safety net (min with it), exactly as the system does.
 * Every candidate is decoded, inverse-filtered and compared with the input.
 *-----------------------------------------------------------------------------*/
static uint8_t *delta_filter(const uint8_t *d, size_t n, int s){
    uint8_t *y = xmalloc(n ? n : 1);
    for(size_t i = 0; i < n; ++i) y[i] = i >= (size_t)s ? (uint8_t)(d[i] - d[i - (size_t)s]) : d[i];
    return y;
}
static int unfilter_ok(const uint8_t *y, size_t n, int s, const uint8_t *d){
    uint8_t *x = xmalloc(n ? n : 1); int ok = 1;
    for(size_t i = 0; i < n; ++i) x[i] = i >= (size_t)s ? (uint8_t)(y[i] + x[i - (size_t)s]) : y[i];
    ok = n == 0 || memcmp(x, d, n) == 0; free(x); return ok;
}
static size_t ent_verified(const uint8_t *d, size_t n, int *fail){
    Buf c; b_init(&c); b_u32(&c, (uint32_t)n); b_push(&c, 2); entropy_compress(d, n, &c);
    Buf r; b_init(&r);
    if(!(decompress(c.p, c.len, &r) && r.len == n && memcmp(r.p, d, n) == 0)) *fail = 1;
    size_t len = c.len; b_free(&c); b_free(&r); return len;
}
static size_t zmin(size_t a, size_t b){ return a < b ? a : b; }
static int filter_compare(int nf, char **files){
    printf("file,n,probe_stride,stored,entropy_raw,lz0,"
           "permatch_auto,filter_lz_auto,permatch_oracle,permatch_oracle_stride,filter_lz_oracle,filter_oracle_stride,"
           "filter_entropy_auto,filter_entropy_oracle,verified\n");
    int rc = 0;
    for(int i = 0; i < nf; ++i){
        size_t n; uint8_t *d = read_all(files[i], &n);
        if(!d){ rc = 2; continue; }
        if(n < 16){ free(d); continue; }
        int fail = 0, P = choose_stride(d, n);
        size_t stored = n + 5, ent = ent_verified(d, n, &fail), lz0 = lz_verified(d, n, 0, &fail);
        size_t pm_auto = lz0, fl_auto = lz0, fe_auto = ent;
        size_t pm_or = lz0, fl_or = lz0, fe_or = ent; int pm_or_s = 0, fl_or_s = 0;
        for(int k = 0; k < NELEM(DELTA_STRIDES); ++k){
            int st = DELTA_STRIDES[k];
            if(!(st > 0 && (size_t)st < (n < 256 ? n : 256))) continue;
            size_t pm = lz_verified(d, n, st, &fail);
            uint8_t *y = delta_filter(d, n, st);
            if(!unfilter_ok(y, n, st, d)) fail = 1;
            size_t fl = lz_verified(y, n, 0, &fail) + 1, fe = ent_verified(y, n, &fail) + 1;
            free(y);
            if(st == P){ pm_auto = zmin(pm, lz0); fl_auto = zmin(fl, lz0); fe_auto = zmin(fe, ent); }
            if(pm < pm_or){ pm_or = pm; pm_or_s = st; }
            if(fl < fl_or){ fl_or = fl; fl_or_s = st; }
            if(fe < fe_or) fe_or = fe;
        }
        const char *b = files[i] + strlen(files[i]);
        while(b > files[i] && b[-1] != '/' && b[-1] != '\\') b--;
        printf("%s,%zu,%d,%zu,%zu,%zu,%zu,%zu,%zu,%d,%zu,%d,%zu,%zu,%s\n", b, n, P, stored, ent, lz0,
               pm_auto, fl_auto, pm_or, pm_or_s, fl_or, fl_or_s, fe_auto, fe_or, fail ? "FAIL" : "ok");
        fflush(stdout);
        if(fail) rc = 4;
        free(d);
    }
    return rc;
}

int main(int argc, char **argv){
    init_cost();
    if(argc < 3 || (argv[1][0] != 't' && argv[1][0] != 'a' && argv[1][0] != 'f' && argc < 4)){
        fprintf(stderr, "usage: %s c <in> <out> | d <in> <out> | t <file> [...] | a <file> [...] | f <file> [...]\n", argv[0]); return 1;
    }
    char mode = argv[1][0];
    if(mode == 'c' || mode == 'd'){
        size_t n; uint8_t *d = read_all(argv[2], &n); if(!d) return 2;
        Buf o; b_init(&o); int ok = 1;
        if(mode == 'c') compress(d, n, &o); else ok = decompress(d, n, &o);
        if(!ok){ fprintf(stderr, "decode failed\n"); return 3; }
        ok = write_all(argv[3], o.p, o.len);
        fprintf(stderr, "%zu -> %zu bytes\n", n, o.len);
        free(d); b_free(&o); return ok ? 0 : 2;
    }
    if(mode == 'a') return ablation(argc - 2, argv + 2);
    if(mode == 'f') return filter_compare(argc - 2, argv + 2);
    if(mode != 't'){ fprintf(stderr, "unknown mode '%c'\n", mode); return 1; }
    int rc = 0;
    for(int i = 2; i < argc; ++i){
        size_t n; uint8_t *d = read_all(argv[i], &n); if(!d){ rc = 2; continue; }
        clock_t t0 = clock(); Buf c; b_init(&c); compress(d, n, &c);
        clock_t t1 = clock(); Buf r; b_init(&r); int ok = decompress(c.p, c.len, &r) && r.len == n && (n == 0 || memcmp(r.p, d, n) == 0);
        clock_t t2 = clock();
        printf("%s: %zu -> %zu (x%.4f) mode=%d enc %.2fs dec %.2fs %s\n", argv[i], n, c.len, c.len ? (double)n / (double)c.len : 0.0,
               c.len > 4 ? c.p[4] : -1, (double)(t1 - t0) / CLOCKS_PER_SEC, (double)(t2 - t1) / CLOCKS_PER_SEC, ok ? "OK" : "LOSSY");
        if(n) printf("  candidates: stored=%zu entropy=%zu", g_info.stored, g_info.entropy);
        if(g_info.has_lz) printf(" lz(stride=%d)=%zu", g_info.stride, g_info.lz);
        if(g_info.has_lz0) printf(" lz(stride=0)=%zu", g_info.lz0);
        printf("  chosen=%s\n", g_info.chosen);
        if(!ok) rc = 4;
        free(d); b_free(&c); b_free(&r);
    }
    return rc;
}
