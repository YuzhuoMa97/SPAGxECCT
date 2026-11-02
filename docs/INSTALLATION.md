# Installation and first run

The [complete official user guide](USER_GUIDE.md) gives checksum verification,
source versus binary directory layouts, step-by-step first runs, all four
methods and 33 runnable feature groups. Start there for a complete tutorial.

R and C++ live side by side in this repository. Installing the historical R
package does not install the independent executable. Windows and compatible Linux
binary users can use a complete Release bundle; source users build the complete
`cpp/` tree.

SPAGxE is a standalone C++17 application. Analysis/tutorials do not call R or
Python; those are optional developer/reference-check dependencies. The
historical R package has a separate installation route.

## Windows precompiled package

Download the complete SPAGxE-1.2.0-windows-x86_64.zip and SHA256SUMS.txt from
[Releases](https://github.com/YuzhuoMa97/SPAGxECCT/releases). These assets exist
only after the maintainer uploads them. Verify SHA256 and extract everything.
In the directory containing run_examples.cmd:

~~~powershell
.\cpp\SPAGxE.exe --version
.\cpp\SPAGxE.exe --help
.\run_examples.cmd
~~~

The prepared binary needs Windows x86-64 with Haswell/AVX2, not R/Python/compiler.
Keep examples and licenses with the EXE. Tutorials write cpp/results/tutorials;
packing refuses existing outputs, so rerun it in a fresh output directory.
Root run_examples.cmd belongs to the binary bundle, not the source checkout.

## Linux precompiled package

Download the complete SPAGxE-1.2.0-linux-x86_64.tar.gz and SHA256SUMS.txt from
[Releases](https://github.com/YuzhuoMa97/SPAGxECCT/releases). The maintainer must
upload these assets before this route is available. In the download directory,
check the Linux archive's own checksum entry, then extract the complete bundle:

~~~sh
grep 'SPAGxE-1.2.0-linux-x86_64.tar.gz$' SHA256SUMS.txt | sha256sum -c -
tar -xzf SPAGxE-1.2.0-linux-x86_64.tar.gz
cd SPAGxE-1.2.0-linux-x86_64
./cpp/SPAGxE --version
./cpp/SPAGxE --help
bash run_examples.sh
~~~

The prepared binary was tested on Ubuntu 24.04.5 under WSL2. It requires an
x86-64 Haswell/AVX2 CPU and compatible dynamic libraries providing at least
GLIBC 2.38, GLIBCXX 3.4.30 and CXXABI 1.3.13. It does not require R, Python
or a compiler; Bash and the compatible system runtime remain necessary.
Older distributions or other CPU targets should use the source-build route
below instead of assuming this binary is portable to every Linux system.
Keep examples and licenses with the executable and run in a writable directory.
Tutorials write cpp/results/tutorials; use a fresh output directory when
repeating a packing example because existing packing outputs are refused.

## Source

Clone or download the repository/tag source ZIP. Read root README_CPP.md,
then build inside cpp/. Keep the complete cpp/ tree, including third_party/,
Makefile, VERSION, examples and tests; src/ alone is insufficient.

Before the first C++ release, select the published `develop` branch explicitly;
the default `main` branch may still contain only the historical R package.
The commands below require the maintainer to have pushed that C++ branch.
After the release is available, replace `--branch develop` with
`--branch cpp-v1.2.0` to obtain the fixed release source. For development
comparisons, record the full tested commit rather than relying on a moving
branch. Code -> Download ZIP likewise uses the selected branch or tag.
See the [development and release procedure](GITHUB_DISTRIBUTION.md).

With C++17 GCC/Clang, a C compiler, GNU Make and Bash:

~~~sh
git clone --branch develop --single-branch https://github.com/YuzhuoMa97/SPAGxECCT.git
cd SPAGxECCT/cpp
mkdir -p tmp results
export TMPDIR="$PWD/tmp"
export TMP="$TMPDIR" TEMP="$TMPDIR"
make -j4
make -j4 test
build/bin/SPAGxE --version
bash examples/tutorials/run_all.sh
~~~

Windows source builds use MSYS2/MinGW or compatible Rtools; the program is
build/bin/SPAGxE.exe. Using an Rtools compiler adds no R runtime dependency.
See [build options](../cpp/README.md). Clean-rebuild after compiler,
architecture/alignment or static-link changes. Bundled sources are compiled
locally without downloads by make. Default x86 builds require Haswell/AVX2.

make install PREFIX=/your/path is optional; it installs the executable in
PREFIX/bin, not all example data. Running build/bin directly is sufficient.
make dist builds a complete bundle in a new dist directory and refuses
overwrite. Custom binaries retain their platform/ISA/runtime requirements.

Linux/macOS routes and CI are supplied; their actual status must be checked
for the released commit. Windows Bash is not a Linux run.

The [tutorial guide](../cpp/examples/tutorials/README.md) has individual
SH/CMD scripts for 01_spagxe, 02_spagxe_plus, 03_spagxemix and
04_spagxemix_local, followed by 05_global_local_cct and 06_local_packing.
Small-example permissive QC is not a real-cohort threshold recommendation.


CI artifacts are development outputs, not a substitute for the complete Release
bundle. The supplied build workflow does not publish a Release automatically.
The maintainer must create the tag/Release and upload the Windows ZIP, Linux
TAR.GZ and SHA256SUMS.txt after verifying the selected source commit.
