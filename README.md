# From Recurrence to Encoded Cost: Selecting Delta Domains for LZ77

# Abstract
Lag differencing can reveal repeated structure that ordinary byte-equality LZ77 cannot detect, but transformed recurrence does not by itself establish a reduction in encoded size. We introduce an evaluation framework that separates five quantities often conflated in transform-selection studies: transformed recurrence, opportunity within a fixed candidate family, selector regret, parser-realizable benefit, and complete encoded cost. We implement this framework in a verified C99 research codec that supports raw and lag-difference LZ77 matching, whole-file delta filtering, multiple entropy coders, candidate-family oracle sweeps, round-trip validation, and targeted ablations. On 62 constructed inputs, a bounded six-byte recurrence selector reduces aggregate mixed-domain archive size from 1,879,861 to 561,096 bytes relative to an entropy selector and finishes only 17 bytes above the fixed-family oracle, demonstrating that it recognizes deliberately embedded transformed structure. This success does not transfer to a prospectively frozen panel of 24 complete files: mixed-domain selection saves only 60 bytes, matches fixed lag one on every file, and requires 1.81 to 1.86 times the raw-policy encoding time. Crucially, the mixed-domain oracle exposes only 271 bytes of aggregate opportunity, showing that the principal limitation is the candidate family rather than selector accuracy alone. Diagnostic ablations identify admission rules, bounded reference visibility, competition with raw matches, and globally coded streams as causes of disagreement between recurrence and archive-size rankings. Whole-file lag filtering yields a distinct 3,855-byte benefit on ten WAV containers through entropy-only coding. The paper contributes a reproducible methodology and verified experimental implementation for determining whether transform selection is limited by recognition, parsing, coding, or the absence of meaningful candidate-family opportunity.

# LZP Research Codec

Research software accompanying the manuscript **“From Recurrence to Encoded Cost: Selecting Delta Domains for LZ77.”**

This repository contains the supplied C99 implementation of a lossless research codec with:

- stored and entropy-only candidates;
- raw-domain and mixed raw/delta-domain LZ77 candidates;
- lag selection by sampled difference entropy;
- context-clustered rANS and adaptive range-coding paths;
- compression and decompression commands;
- round-trip verification;
- delta-domain ablations; and
- per-match versus whole-file delta-filter comparisons.

## Repository layout

```text
.
├── src/lzp_codec.c          Supplied C99 codec and experimental modes
├── scripts/smoke_test.sh    Build and round-trip smoke test
├── tests/data/sample.txt    Small test input generated for this repository
├── docs/REPRODUCING.md      Reproduction and artifact guidance
├── docs/MANUSCRIPT.md       Suggested paper citation and availability text
├── CITATION.cff             GitHub citation metadata
├── LICENSE                  CC0-1.0 dedication
├── Makefile                 Portable GCC/Clang build targets
└── SHA256SUMS               Integrity hashes for tracked research inputs
```

## Build

A C99 compiler and the system math library are required.

```bash
make
```

Equivalent command:

```bash
gcc -O2 -std=c99 -Wall -Wextra -pedantic -o lzp_codec src/lzp_codec.c -lm
```

The strict C99 mode is intentional because the supplied implementation documents floating-point behavior needed for compatibility with the original Python implementation.

## Usage

```bash
# Compress
./lzp_codec c input.bin output.lzp

# Decompress
./lzp_codec d output.lzp restored.bin

# Compress, decompress, verify, and print candidate information
./lzp_codec t input.bin [more-files ...]

# Delta-domain ablation report as CSV
./lzp_codec a input.bin [more-files ...]

# Per-match delta versus whole-file delta-filter report as CSV
./lzp_codec f input.bin [more-files ...]
```

## License

The supplied source header designates the implementation as “Public domain / CC0.” This repository therefore includes the CC0 1.0 Universal dedication. See [`LICENSE`](LICENSE).

## Research-use disclaimer

This is an experimental research codec, not a production compression library. Validate archives independently and retain original data.
