# Contributing

Contributions that improve correctness, reproducibility, portability, documentation, or test coverage are welcome.

## Before opening a pull request

1. Build with `make clean && make`.
2. Run `make test`.
3. Preserve the archive format unless the change explicitly introduces a documented format version.
4. Add or update tests for behavior changes.
5. Document any change that can alter candidate ordering, tie-breaking, floating-point accumulation, or compressed bytes.

## Reproducibility-sensitive code

The implementation intentionally preserves specific sorting, tie-breaking, hashing, and numerical behavior. Refactoring these areas can change compressed output even when decompression remains correct. Pull requests should report whether generated archives remain byte-identical on the regression corpus.

## Reporting issues

Include the compiler and version, operating system, exact command, input hash when shareable, observed output, and expected output. Do not attach restricted or confidential datasets.
