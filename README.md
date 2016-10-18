To do:
  - Integrate all code - c-line parser, samparser.cpp, bwa call
    - Integrated, not tested
  - Build suffix array reader - pre-gen input for bio construction
  - Remove experimental code
  - Remove hard code
  - Simplify code after suffix array construction


Optimizations (biological):
  - Streamline biological performance workflow (python analysis)
  - Reduction of false positives (consensus sequence manipulation, construction)
  - Increase true positives

Optimizations (resources):
  - 2-bit compression scheme
  - Work out way to recreate div-suf-sort with 256 byte alphabet
  - Initial scan, removing non-mutated sequences from cancer data set by
    alignment
  - either string graph representation, or alignment compression
