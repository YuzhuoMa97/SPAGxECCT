# Incorporated source and dependency notices

SPAGxE algorithms and new standalone changes: Yuzhuo Ma, yuzhuoma@amss.ac.cn.
The method's published authorship is retained in CITATION_CPP.cff.

New SPAGxE work includes the paper-defined interaction policy, complete-CGF
workflow, AF/input/numerical corrections, local G/H analysis and packing,
standalone interfaces, examples and validation. Historical R files in this
repository keep their existing contents and license scope.

| Component | Location | License |
|---|---|---|
| Incorporated GPL application code | cpp/src, Makefile, tests | GPL-3.0-or-later; cpp/LICENSE |
| Eigen 5.0.0 | cpp/third_party/eigen-5.0.0 | MPL-2.0 and included notices |
| Boost 1.80.0 | cpp/third_party/boost-1.80.0 | Boost Software License 1.0 |
| pgenlib / libdeflate | cpp/third_party/plink2-a.6.33 | Included LGPL-3-or-later / MIT notices |
| BGEN library | cpp/third_party/bgen-1.2.0 | Boost Software License 1.0 |
| htslib | cpp/third_party/htslib-1.23.1 | MIT and included notices |
| zlib | cpp/third_party/zlib-1.3.2 | zlib license |
| zstd | cpp/third_party/zstd-1.5.7 | Included BSD/GPL choices |
| SQLite | cpp/third_party/sqlite-3.53.2 | Public domain dedication |

The public examples/1kg dataset is retained from the same upstream snapshot.
Other bundled association examples are synthetic SPAGxE data.
No private external-project source, manuscripts or real study datasets are
included in the new standalone distribution.

Windows source builds additionally use POSIX-regex libraries supplied by the
selected MinGW toolchain. A custom source build may have runtime requirements
different from the tested, statically linked Windows release bundle.
