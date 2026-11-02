# Standalone SPAGxE C++ 1.2.0

**Complete official teaching manual:** [download, installation, every current workflow, outputs and troubleshooting](docs/USER_GUIDE.md).

Run the extended feature tutorial after the six introductory examples: `examples/tutorials/07_feature_walkthrough.sh` or `examples/tutorials/07_feature_walkthrough.cmd` from `cpp/`. See the guide for executable selection and fresh output directories.

SPAGxEmixPlus is planned for a future release and is not included in the current C++ distribution.
For ancestry-specific files, start with [local input preparation](docs/LOCAL_INPUT_PREPARATION.md).
Development changes are tested on `develop` before merge to `main`; follow
the [development and release workflow](docs/GITHUB_DISTRIBUTION.md).

Maintainer: **Yuzhuo Ma, yuzhuoma@amss.ac.cn**.
Four methods: SPAGxE, SPAGxE+, SPAGxEmix and SPAGxEmix_local.

## Official C++ distribution and support

This repository is the sole official standalone C++ distribution of the
SPAGxE method family designated and maintained by **Yuzhuo Ma**, the method
designer and author of the original SPAGxE software.
For questions, bug reports, installation help or result interpretation,
contact **[yuzhuoma@amss.ac.cn](mailto:yuzhuoma@amss.ac.cn)**.

The historical R package remains available for reproducibility and compatibility.
Published-paper coauthorship and incorporated-source authorship and licenses
remain as recorded in CITATION_CPP.cff, OWNERSHIP.md and THIRD_PARTY_NOTICES.md.

Read [installation](docs/INSTALLATION.md), [the full CLI/build manual](cpp/README.md),
and [annotated method tutorials](cpp/examples/tutorials/README.md).
The analysis executable does not call R or Python.

Source users clone/download the complete repository and build inside cpp/:

~~~sh
cd cpp
mkdir -p tmp results
export TMPDIR="$PWD/tmp"
export TMP="$TMPDIR" TEMP="$TMPDIR"
make -j4
make -j4 test
build/bin/SPAGxE --version
bash examples/tutorials/run_all.sh
~~~

Windows users may instead download the complete SPAGxE-1.2.0-windows-x86_64.zip
asset from [Releases](https://github.com/YuzhuoMa97/SPAGxECCT/releases)
after it is uploaded. Extract it, run .\cpp\SPAGxE.exe --version and
.\run_examples.cmd. The prepared binary requires Windows x86-64 and
Haswell/AVX2; it needs no R, Python or compiler at runtime.
Windows source builds use MinGW/MSYS2/Rtools.

There is one statistical policy: paper-defined weights and guarded partial-normal
CGF acceleration for all four methods, with automatic full-CGF protection.
This adds an approximation; it does not promise unchanged historical P values.
Read [partial SPA](docs/PARTIAL_SPA.md),
[method policy and limits](docs/METHOD_POLICY.md) and
[validation](docs/validation/README.md) before interpreting results.
Windows has been executed locally; check actual Actions for other platforms.

The historical R package remains version 1.1.0 with its original dependencies.
R installation does not install the C++ executable. Existing R algorithms,
DESCRIPTION, NAMESPACE and root LICENSE are preserved.
The C++ application uses cpp/LICENSE and retained
[component notices](THIRD_PARTY_NOTICES.md).
