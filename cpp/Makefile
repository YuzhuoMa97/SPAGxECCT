# ── Toolchain ─────────────────────────────────────────────────────────────────
# External build environments and CI may override the compiler defaults.
# The supplied workflow builds and tests; Release assets are uploaded separately.
ifeq ($(origin CXX),default)
CXX = g++
endif
CXX ?= g++
ifeq ($(origin CC),default)
CC = gcc
endif
CC ?= gcc

# ── Version ───────────────────────────────────────────────────────────────────
# VERSION in this cpp/ directory supplies the SPAGXE_VERSION compile macro.
# Keep the version, documentation and maintainer-created Release assets aligned.
# The first standalone C++ release uses tag cpp-v1.2.0 in the parent repository.
# CI build artifacts are not complete Release bundles and are not auto-published.
SPAGXE_VERSION ?= $(strip $(shell cat VERSION 2>/dev/null || echo 0.0.0+unknown))

# ── Platform detection ────────────────────────────────────────────────────────
# UNAME_S is empty on plain Windows cmd; MSYS2/MinGW sets it to MINGW* or MSYS*.
UNAME_S := $(shell uname -s 2>/dev/null)

ifeq ($(filter MINGW% MSYS%,$(UNAME_S)),)
  ifeq ($(filter Darwin%,$(UNAME_S)),)
    # ── Linux ─────────────────────────────────────────────────────────────────
    PLATFORM  := linux
    EXE       :=
    PLATFORM_LDLIBS :=
    # -pipe: pass data between gcc stages via pipes (no temp files).
    PLATFORM_FLAGS  := -pipe
    SHELL     := /bin/bash
  else
    # ── macOS ─────────────────────────────────────────────────────────────────
    PLATFORM  := macos
    EXE       :=
    PLATFORM_LDLIBS :=
    PLATFORM_FLAGS  := -pipe
    SHELL     := /bin/bash
  endif
else
  # ── Windows (MSYS2 / MinGW / Rtools) ─────────────────────────────────────
  PLATFORM  := windows
  EXE       := .exe
  # ws2_32: Winsock2 (htslib hfile.c socket I/O)
  # shell32: CommandLineToArgvW (UTF-8-preserving native Windows CLI paths).
  # regex:  POSIX regex (htslib hts_expr.c).  Static-link chain on
  # MSYS2/MinGW:
  #   libregex.a → tre_regcomp / tre_regfree / tre_regerror / tre_regexec
  #   libtre.a   → libintl_gettext (i18n of tre's error strings)
  #   libintl.a  → libiconv_open / libiconv_close / libiconv
  # In dynamic-mode builds the libregex DLL has the entire chain
  # baked in, so the explicit `-ltre -lintl -liconv` block is a
  # no-op; in static-mode builds each .a is required.  Wrapping the
  # block in -Bstatic/-Bdynamic forces these three to be picked up
  # from .a archives regardless of the surrounding global link mode
  # so the resulting .exe is fully self-contained.
  PLATFORM_LDLIBS := -lws2_32 -lshell32 -lregex -Wl,-Bstatic -ltre -lintl -liconv
  PLATFORM_FLAGS  := -pipe
  SHELL     := /usr/bin/bash
  # Redirect linker temp dir to project-local tmp/ (Windows path format).
  WIN_TMP   := $(shell cygpath -w $(CURDIR)/tmp)
  export TMP  := $(WIN_TMP)
  export TEMP := $(WIN_TMP)
endif

# ── Architecture detection ────────────────────────────────────────────────────
# IS_X86 = 1 on x86_64; otherwise 0 (covers arm64/aarch64 macOS M-series, etc.).
UNAME_M := $(shell uname -m 2>/dev/null)
ifeq ($(UNAME_M),x86_64)
  IS_X86 := 1
else
  IS_X86 := 0
endif

# ── AVX2 detection (x86 only) ─────────────────────────────────────────────────
# Apple clang on ARM warns "-mavx2 unused" instead of erroring, so the probe
# would falsely succeed.  Restrict the probe to x86 outright.
ifeq ($(IS_X86),1)
  AVX2_OK := $(shell TMP_AVX=$$(mktemp -u 2>/dev/null || echo __avx2_test__); \
    echo 'int main(){return 0;}' | $(CXX) -x c++ -mavx2 - -o "$$TMP_AVX" 2>/dev/null \
    && { rm -f "$$TMP_AVX"; echo 1; } || rm -f "$$TMP_AVX")
  ifeq ($(AVX2_OK),1)
    SIMD_FLAGS := -mavx2 -mbmi -mbmi2 -mlzcnt -mfma
  else
    SIMD_FLAGS :=
  endif
else
  SIMD_FLAGS :=
endif

# Architecture baseline for application sources.
# Selected kernels dispatch to target-specific variants at runtime. This does
# not make the whole executable compatible with CPUs below its build baseline:
# application auto-vectorization and third-party code may require AVX2.
ifeq ($(IS_X86),1)
  # Default x86 baseline: Haswell (AVX2/FMA/BMI2), matching the prepared binary.
  # A host-specific build may override SPAGXE_MARCH, followed by a clean rebuild.
  # Lowering this flag alone does not lower the separate third-party SIMD_FLAGS.
  # Target-specific AVX-512 kernels may still be selected on capable hosts.
  SPAGXE_MARCH ?= -march=haswell
else
  # ARM (arm64 / aarch64): rely on default tuning + scalar fallbacks.
  SPAGXE_MARCH ?=
endif

# ── Compiler & linker flags ───────────────────────────────────────────────────
# Naming convention:
#   TP_CXXFLAGS / TP_CFLAGS  : internal defaults for third-party (zlib, zstd,
#                              libdeflate, pgenlib, bgen, htslib).  Include
#                              SIMD_FLAGS so the C/C++ deps get AVX2 codegen.
#   SPAGXE_CXXFLAGS         : application defaults including SPAGXE_MARCH;
#                            selected kernels also use runtime SIMD dispatch.
#   LINK_FLAGS               : internal defaults for the final link step.
#
# The standard $(CPPFLAGS) / $(CXXFLAGS) / $(CFLAGS) / $(LDFLAGS) names are
# left for the build environment (external sysroots, hardening flags such
# as -fstack-protector-strong / -D_FORTIFY_SOURCE=2 / -fdebug-prefix-map=…,
# and rpath/link-search additions) to fill in.  They are appended after the
# internal flag sets so external overrides win for any conflicting option.
TP_CXXFLAGS := -std=c++17 -O3 -DNDEBUG $(PLATFORM_FLAGS) \
            -ffunction-sections -fdata-sections \
            -funroll-loops $(SIMD_FLAGS) \
            -Wall -Wextra -Wno-unused-parameter -Wno-sign-compare
SPAGXE_CXXFLAGS := -std=c++17 -O3 -DNDEBUG -DEIGEN_MAX_ALIGN_BYTES=16 -DEIGEN_MAX_STATIC_ALIGN_BYTES=16 $(SPAGXE_MARCH) $(PLATFORM_FLAGS) \
            -ffunction-sections -fdata-sections \
            -funroll-loops \
            -DSPAGXE_VERSION=\"$(SPAGXE_VERSION)\" \
            -Wall -Wextra -Wno-unused-parameter -Wno-sign-compare \
            -Wno-maybe-uninitialized
TP_CFLAGS   := -O3 -DNDEBUG $(PLATFORM_FLAGS) -ffunction-sections -fdata-sections $(SIMD_FLAGS)
# libstdc++/libgcc linkage: dynamic by default (works everywhere; requires the
# system's libstdc++.so.6/libgcc_s.so.1 at runtime — both are part of any
# Linux/Windows g++ install).  Override with STATIC_LIBS="-static-libstdc++
# -static-libgcc" if you have the libstdc++-static archive installed and want
# a binary portable to older systems.  Apple clang on macOS doesn't accept
# these flags, so the override is Linux/Windows-only.
STATIC_LIBS ?=
ifeq ($(PLATFORM),windows)
  # Restore dynamic lookup only for a dynamic build. GCC translates -static
  # into an early linker switch; a later -Bdynamic would override it for
  # implicit C++/thread runtime libraries and leave unexpected DLL imports.
  ifeq ($(filter -static,$(STATIC_LIBS) $(LDFLAGS) $(CXXFLAGS)),)
    PLATFORM_LDLIBS += -Wl,-Bdynamic
  endif
endif
ifeq ($(PLATFORM),macos)
  LINK_FLAGS := -Wl,-dead_strip -lpthread $(PLATFORM_LDLIBS) $(STATIC_LIBS)
else
  # -lstdc++fs is required by libstdc++ before GCC 9 to resolve
  # std::filesystem symbols. Newer toolchains expose an empty stub library
  # so the flag is harmless there.
  LINK_FLAGS := -Wl,--gc-sections -lpthread -lstdc++fs $(PLATFORM_LDLIBS) $(STATIC_LIBS)
endif

# ── Directories ───────────────────────────────────────────────────────────────
SRC_DIR      := src
BUILD_DIR    := build
ZLIB_DIR     := third_party/zlib-1.3.2
ZSTD_DIR     := third_party/zstd-1.5.7/lib
PGENLIB_DIR  := third_party/plink2-a.6.33
DEFLATE_DIR  := $(PGENLIB_DIR)/libdeflate
BGEN_DIR     := third_party/bgen-1.2.0
HTSLIB_DIR   := third_party/htslib-1.23.1
SQLITE_DIR   := third_party/sqlite-3.53.2

# ── Output binary ─────────────────────────────────────────────────────────────
BIN := $(BUILD_DIR)/bin/SPAGxE$(EXE)

# ── Include paths ─────────────────────────────────────────────────────────────
# Boost.Math and Eigen are header-only — no -l flags required.
INCLUDES := \
    -I$(SRC_DIR) \
    -Ithird_party/eigen-5.0.0 \
    -Ithird_party/boost-1.80.0 \
    -I$(ZLIB_DIR) \
    -I$(ZSTD_DIR) \
    -I$(PGENLIB_DIR)/include \
    -I$(BGEN_DIR)/genfile/include \
    -I$(HTSLIB_DIR) \
    -I$(HTSLIB_DIR)/htscodecs \
    -I$(SQLITE_DIR) \
    -I$(DEFLATE_DIR)

# ── Source discovery ──────────────────────────────────────────────────────────
# Pure-Make recursive wildcard: avoids calling the shell 'find' command,
# which can resolve to Windows FIND.EXE instead of the Unix utility.
rwildcard = $(foreach d,$(wildcard $(addsuffix /*,$(1))),\
                $(call rwildcard,$(d),$(2)) $(filter $(subst *,%,$(2)),$(d)))

# ── Application sources (C++) ────────────────────────────────────────────
SRCS := $(SRC_DIR)/main.cpp $(wildcard $(SRC_DIR)/cli/*.cpp) $(SRC_DIR)/engine/marker.cpp $(wildcard $(SRC_DIR)/io/*.cpp) $(wildcard $(SRC_DIR)/genotype/*.cpp) $(wildcard $(SRC_DIR)/interaction/*.cpp) $(SRC_DIR)/allele_model/frequency_model.cpp $(wildcard $(SRC_DIR)/local/*.cpp) $(wildcard $(SRC_DIR)/util/*.cpp)
OBJS := $(patsubst $(SRC_DIR)/%.cpp, $(BUILD_DIR)/%.o, $(SRCS))
DEPS := $(OBJS:.o=.d)

# ── zlib (C) ──────────────────────────────────────────────────────────────────
ZLIB_SRCS := $(wildcard $(ZLIB_DIR)/*.c)
ZLIB_OBJS := $(patsubst $(ZLIB_DIR)/%.c, $(BUILD_DIR)/zlib/%.o, $(ZLIB_SRCS))

# ── zstd (C) — shared by pgenlib and bgen ─────────────────────────────────────
ZSTD_SRCS := $(wildcard $(ZSTD_DIR)/common/*.c) \
             $(wildcard $(ZSTD_DIR)/compress/*.c) \
             $(wildcard $(ZSTD_DIR)/decompress/*.c)
ZSTD_OBJS := $(patsubst $(ZSTD_DIR)/%.c, $(BUILD_DIR)/zstd/%.o, $(ZSTD_SRCS))

# ── libdeflate (C) — used by pgenlib ──────────────────────────────────────────
DEFLATE_SRCS := $(wildcard $(DEFLATE_DIR)/lib/*.c)
# Per-architecture CPU-feature detection.  Each *_cpu_features.c is guarded
# internally by arch #ifdefs, so compiling both on any target is harmless (the
# non-matching one becomes an empty translation unit); the matching one supplies
# the runtime-dispatch symbols (libdeflate_{x86,arm}_cpu_features) referenced by
# adler32.c / crc32.c.  Omitting the ARM file breaks the aarch64 link.
ifneq ($(wildcard $(DEFLATE_DIR)/lib/x86/x86_cpu_features.c),)
  DEFLATE_SRCS += $(DEFLATE_DIR)/lib/x86/x86_cpu_features.c
endif
ifneq ($(wildcard $(DEFLATE_DIR)/lib/arm/arm_cpu_features.c),)
  DEFLATE_SRCS += $(DEFLATE_DIR)/lib/arm/arm_cpu_features.c
endif
DEFLATE_OBJS := $(patsubst $(DEFLATE_DIR)/lib/%.c, $(BUILD_DIR)/libdeflate/%.o, $(DEFLATE_SRCS))

# ── pgenlib (C++ .cc files + SFMT.c) ─────────────────────────────────────────
PGEN_CC_SRCS := $(wildcard $(PGENLIB_DIR)/include/*.cc)
PGEN_CC_OBJS := $(patsubst $(PGENLIB_DIR)/include/%.cc, $(BUILD_DIR)/pgenlib/%.o, $(PGEN_CC_SRCS))
PGEN_C_SRCS  := $(PGENLIB_DIR)/include/SFMT.c
PGEN_C_OBJS  := $(BUILD_DIR)/pgenlib/SFMT.o
PGEN_OBJS    := $(PGEN_CC_OBJS) $(PGEN_C_OBJS)

# ── bgen (C++) — exclude View.cpp (needs boost::filesystem + IndexQuery) ─────
BGEN_SRCS := $(BGEN_DIR)/src/bgen.cpp \
             $(BGEN_DIR)/src/zlib.cpp \
             $(BGEN_DIR)/src/MissingValue.cpp
BGEN_OBJS := $(patsubst $(BGEN_DIR)/src/%.cpp, $(BUILD_DIR)/bgen/%.o, $(BGEN_SRCS))

# ── sqlite (C) — read-only .bgi (bgenix) index access for BGEN ─────────────────
SQLITE_SRCS := $(SQLITE_DIR)/sqlite3.c
SQLITE_OBJS := $(BUILD_DIR)/sqlite/sqlite3.o

# ── htslib (C) ────────────────────────────────────────────────────────────────
HTSLIB_SRCS    := $(wildcard $(HTSLIB_DIR)/*.c)
HTSCODEC_SRCS  := $(wildcard $(HTSLIB_DIR)/htscodecs/htscodecs/*.c)
# os/rand.c is #included by hts_os.c on Windows — not compiled separately.
HTSLIB_OBJS    := $(patsubst $(HTSLIB_DIR)/%.c, $(BUILD_DIR)/htslib/%.o, $(HTSLIB_SRCS))
HTSCODEC_OBJS  := $(patsubst $(HTSLIB_DIR)/htscodecs/htscodecs/%.c, \
                     $(BUILD_DIR)/htslib/htscodecs/%.o, $(HTSCODEC_SRCS))

# ── All objects ───────────────────────────────────────────────────────────────
ALL_OBJS := $(OBJS) $(ZLIB_OBJS) $(ZSTD_OBJS) $(DEFLATE_OBJS) \
            $(PGEN_OBJS) $(BGEN_OBJS) \
            $(HTSLIB_OBJS) $(HTSCODEC_OBJS) $(SQLITE_OBJS)

# ── Default target ────────────────────────────────────────────────────────────
.PHONY: all clean run

all: $(BIN)

# ── Link ──────────────────────────────────────────────────────────────────────
$(BIN): $(ALL_OBJS) | $(BUILD_DIR) tmp
	@mkdir -p $(@D)
	$(CXX) $(TP_CXXFLAGS) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) $^ -o $@ $(LINK_FLAGS) $(LDLIBS)

# ── Directory creation ────────────────────────────────────────────────────────
$(BUILD_DIR) tmp:
	mkdir -p $@

# ── Compile application (C++) ────────────────────────────────────────────────────────
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp | tmp
	@mkdir -p $(@D)
	$(CXX) $(SPAGXE_CXXFLAGS) $(CPPFLAGS) $(CXXFLAGS) -MMD -MP $(INCLUDES) -c $< -o $@

-include $(DEPS)

# ── Compile zlib (C) ──────────────────────────────────────────────────────────
$(BUILD_DIR)/zlib/%.o: $(ZLIB_DIR)/%.c | tmp
	@mkdir -p $(@D)
	$(CC) $(TP_CFLAGS) $(CPPFLAGS) $(CFLAGS) -I$(ZLIB_DIR) -c $< -o $@

# ── Compile zstd (C) ──────────────────────────────────────────────────────────
$(BUILD_DIR)/zstd/%.o: $(ZSTD_DIR)/%.c | tmp
	@mkdir -p $(@D)
	$(CC) $(TP_CFLAGS) $(CPPFLAGS) $(CFLAGS) -DZSTD_DISABLE_ASM -I$(ZSTD_DIR) -I$(ZSTD_DIR)/common -c $< -o $@

# ── Compile libdeflate (C) ────────────────────────────────────────────────────
$(BUILD_DIR)/libdeflate/%.o: $(DEFLATE_DIR)/lib/%.c | tmp
	@mkdir -p $(@D)
	$(CC) $(TP_CFLAGS) $(CPPFLAGS) $(CFLAGS) -I$(DEFLATE_DIR) -c $< -o $@

# ── Compile pgenlib (.cc → C++) ───────────────────────────────────────────────
#  IGNORE_BUNDLED_ZSTD: use our vendored zstd via -I instead of ../zstd/ path.
PGEN_CXXFLAGS := $(TP_CXXFLAGS) -DIGNORE_BUNDLED_ZSTD \
    -Wno-sign-compare -Wno-unused-function -Wno-missing-field-initializers \
    -Wno-maybe-uninitialized

$(BUILD_DIR)/pgenlib/%.o: $(PGENLIB_DIR)/include/%.cc | tmp
	@mkdir -p $(@D)
	$(CXX) $(PGEN_CXXFLAGS) $(CPPFLAGS) $(CXXFLAGS) -I$(PGENLIB_DIR)/include -I$(ZSTD_DIR) \
	    -I$(DEFLATE_DIR) -c $< -o $@

# SFMT.c is plain C
$(BUILD_DIR)/pgenlib/SFMT.o: $(PGENLIB_DIR)/include/SFMT.c | tmp
	@mkdir -p $(@D)
	$(CC) $(TP_CFLAGS) $(CPPFLAGS) $(CFLAGS) -I$(PGENLIB_DIR)/include -c $< -o $@

# ── Compile bgen (C++) ────────────────────────────────────────────────────────
BGEN_CXXFLAGS := $(TP_CXXFLAGS) -Wno-sign-compare -Wno-unused-variable

$(BUILD_DIR)/bgen/%.o: $(BGEN_DIR)/src/%.cpp | tmp
	@mkdir -p $(@D)
	$(CXX) $(BGEN_CXXFLAGS) $(CPPFLAGS) $(CXXFLAGS) -I$(BGEN_DIR)/genfile/include \
	    -I$(ZSTD_DIR) -I$(ZLIB_DIR) -c $< -o $@

# ── Compile sqlite (C) — read-only .bgi index; lean, no extension loading ─────
# SQLITE_OMIT_LOAD_EXTENSION avoids a libdl dependency; THREADSAFE=1 keeps the
# one-shot index read (main thread, BgenCatalog construction) safe in the otherwise
# multi-threaded process.  No -march=native (uses third-party SIMD_FLAGS only).
# -Wno-stringop-overread silences a GCC false positive in sqlite3ColumnSetColl
# (the inlined sqlite3Strlen30 -> strlen call): GCC cannot prove the collation-
# name pointer is non-null at that site and warns about a "region of size 0".
SQLITE_CFLAGS := $(TP_CFLAGS) \
    -DSQLITE_OMIT_LOAD_EXTENSION -DSQLITE_THREADSAFE=1 \
    -DSQLITE_DEFAULT_MEMSTATUS=0 -DSQLITE_OMIT_DEPRECATED \
    -Wno-stringop-overread

$(BUILD_DIR)/sqlite/%.o: $(SQLITE_DIR)/%.c | tmp
	@mkdir -p $(@D)
	$(CC) $(SQLITE_CFLAGS) $(CPPFLAGS) $(CFLAGS) -I$(SQLITE_DIR) -c $< -o $@

# ── Compile htslib (C) ────────────────────────────────────────────────────────
HTSLIB_CFLAGS := $(TP_CFLAGS) -DHAVE_CONFIG_H \
    -I$(HTSLIB_DIR) -I$(HTSLIB_DIR)/htscodecs -I$(ZLIB_DIR) \
    -Wno-sign-compare

$(BUILD_DIR)/htslib/%.o: $(HTSLIB_DIR)/%.c | tmp
	@mkdir -p $(@D)
	$(CC) $(HTSLIB_CFLAGS) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(BUILD_DIR)/htslib/htscodecs/%.o: $(HTSLIB_DIR)/htscodecs/htscodecs/%.c | tmp
	@mkdir -p $(@D)
	$(CC) $(HTSLIB_CFLAGS) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

# ── Install ───────────────────────────────────────────────────────────────────
# Standard PREFIX / DESTDIR semantics for local or CI installation
# and ad-hoc local installations (`make install PREFIX=$HOME/.local`) can
# stage the binary into a chroot via DESTDIR while keeping the runtime path
# under PREFIX.  Both `install -d` and `install -m` are supported by GNU
# coreutils and BSD `install` (macOS), avoiding the GNU-only `-D` shortcut.
PREFIX  ?= /usr/local
DESTDIR ?=

.PHONY: install
install: $(BIN)
	install -d $(DESTDIR)$(PREFIX)/bin
	install -m 0755 $(BIN) $(DESTDIR)$(PREFIX)/bin/SPAGxE$(EXE)

# ── Dist (release packaging) ──────────────────────────────────────────────────
# A complete runnable bundle, separate from the source checkout:
#   cpp/SPAGxE[.exe], cpp/examples, cpp/README.md, cpp/LICENSE, cpp/VERSION,
#   docs, dependency notices, and root run_examples.sh / run_examples.cmd.
# This target does not invoke a compiler other than the normal $(BIN) prerequisite,
# and the extracted examples need no R, Python or compiler.  Existing output
# directories/archives are refused; choose a new DIST_DIR for another bundle.
# No release CI workflow is assumed.  CI build jobs do not by themselves prove
# that a platform-specific distribution has passed the extracted-example checks.
# The exact input binary is copied by default. DIST_STRIP=1 opts into strip;
# recheck the resulting binary hash if stripping a previously validated binary.
STRIP   ?= strip
DIST_STRIP ?= 0
DIST_DIR ?= dist
DIST_PKG := SPAGxE-$(SPAGXE_VERSION)-$(PLATFORM)-$(UNAME_M)
# Pure-Make discovery avoids a Windows FIND.EXE collision in a native PATH.
DIST_LICENSE_PATTERNS := LICENSE% COPYING% COPYRIGHT% NOTICE% AUTHORS% DISCLAIMER% license% copying% copyright% notice% authors% disclaimer%
DIST_LICENSE_FILES = $(foreach f,$(call rwildcard,third_party,*),$(if $(filter $(DIST_LICENSE_PATTERNS),$(notdir $(f))),$(f)))

.PHONY: dist
dist: $(BIN)
	@set -eu; \
	  case "$(DIST_PKG)" in ""|"."|".."|*[!A-Za-z0-9._-]*) echo "Unsafe distribution package name" >&2; exit 1;; esac; \
	  bundle="$(DIST_DIR)/$(DIST_PKG)"; archive="$(DIST_DIR)/$(DIST_PKG).tar.gz"; \
	  if [ -e "$$bundle" ] || [ -L "$$bundle" ] || [ -e "$$archive" ] || [ -L "$$archive" ]; then \
	    echo "Refusing to overwrite an existing distribution: $$bundle or $$archive" >&2; exit 1; \
	  fi; \
	  case "$(DIST_STRIP)" in 0|1) ;; *) echo "DIST_STRIP must be 0 or 1" >&2; exit 1;; esac; \
	  if [ ! -d ../docs ]; then echo "Distribution requires the matching repository docs directory" >&2; exit 1; fi; \
	  mkdir -p "$(DIST_DIR)"; mkdir "$$bundle"; mkdir -p "$$bundle/cpp" "$$bundle/licenses/third_party"; \
	  cp "$(BIN)" "$$bundle/cpp/SPAGxE$(EXE)"; \
	  if [ "$(DIST_STRIP)" = 1 ]; then $(STRIP) "$$bundle/cpp/SPAGxE$(EXE)"; fi; \
	  cp README.md LICENSE VERSION "$$bundle/cpp/"; cp LICENSE "$$bundle/LICENSE"; \
	  cp -R examples "$$bundle/cpp/examples"; \
	  cp -R ../docs "$$bundle/docs"; \
	  for name in OWNERSHIP.md THIRD_PARTY_NOTICES.md CITATION.cff; do \
	    if [ -f "$$name" ]; then cp "$$name" "$$bundle/$$name"; \
	    elif [ -f "../$$name" ]; then cp "../$$name" "$$bundle/$$name"; fi; \
	  done; \
	  for name in README_CPP.md CITATION_CPP.cff CHANGELOG_CPP.md; do \
	    if [ -f "../$$name" ]; then cp "../$$name" "$$bundle/$$name"; fi; \
	  done; \
	  for file in $(foreach f,$(DIST_LICENSE_FILES),"$(f)"); do \
	    if [ -f "$$file" ]; then rel="$${file#third_party/}"; target="$$bundle/licenses/third_party/$$rel"; \
	      mkdir -p "$${target%/*}"; cp "$$file" "$$target"; fi; \
	  done; \
	  cp examples/run_release_examples.sh "$$bundle/run_examples.sh"; \
	  cp examples/run_release_examples.cmd "$$bundle/run_examples.cmd"; \
	  chmod +x "$$bundle/run_examples.sh"; \
	  printf '%s\n' '# SPAGxE $(SPAGXE_VERSION) binary bundle' '' \
	    'Extract the complete folder to a writable working drive.' \
	    'Run ./run_examples.sh on POSIX or .\run_examples.cmd in Windows cmd/PowerShell.' \
	    'Examples write cpp/results/tutorials inside this extracted folder.' \
	    'Analysis and examples need the executable and shell; no R, Python or compiler is required.' \
	    'Use a binary for your operating system and architecture. Default x86 builds require Haswell/AVX2.' \
	    'Any native runtime libraries required by a custom build remain platform requirements.' \
	    'Source-oriented links in cpp/README.md and docs refer to the matching source checkout.' \
	    'The bundle root LICENSE is the C++ GPL; the existing R repository license is not changed.' \
	    'Maintainer: Yuzhuo Ma, yuzhuoma@amss.ac.cn. See included attribution and dependency notices.' \
	    > "$$bundle/README_RELEASE.md"; \
	  tar -czf "$$archive" -C "$(DIST_DIR)" "$(DIST_PKG)"; \
	  echo "Created $$archive"

# ── Tests and benchmarks ──────────────────────────────────────────────────────
# Developer-only targets.  `tests/` is deliberately outside the rwildcard over
# $(SRC_DIR), so nothing here is linked into the shipped binary and the
# "binary is the deliverable" property is preserved.
#
# Each test/benchmark is a self-contained translation unit with its own main().
# Those that exercise only header-level code (the SPA kernels, the reference
# CGF) link against nothing; those that need compiled application objects list them
# explicitly in TEST_<name>_OBJS below.  There is no test-runner library and no
# third-party framework: tests/tinytest.hpp is ~150 lines of header.
#
#   make test    build and run every test, stop at the first failing binary
#   make bench   build and run every benchmark
#
# Tests are compiled with the same SPAGXE_CXXFLAGS as src/ so that the SIMD
# dispatch, -march setting and optimization level under test match the shipped
# configuration.  Assertions are wanted, so -DNDEBUG is filtered out.
TEST_DIR := tests
TEST_BUILD_DIR := $(BUILD_DIR)/tests
TEST_SRCS := $(wildcard $(TEST_DIR)/*_test.cpp)
TEST_BINS := $(patsubst $(TEST_DIR)/%.cpp,$(TEST_BUILD_DIR)/%$(EXE),$(TEST_SRCS))
TEST_CXXFLAGS := $(filter-out -DNDEBUG,$(SPAGXE_CXXFLAGS))
TESTOBJS_spa_cgf_test := $(BUILD_DIR)/util/spa_cgf.o
TESTOBJS_binomial_cgf_test := $(BUILD_DIR)/util/spa_cgf.o
TESTOBJS_log10p_test := $(BUILD_DIR)/util/probability_math.o $(BUILD_DIR)/util/wald.o $(BUILD_DIR)/genotype/equilibrium.o $(BUILD_DIR)/pgenlib/plink2_stats.o $(BUILD_DIR)/pgenlib/plink2_base.o $(BUILD_DIR)/pgenlib/SFMT.o
TESTOBJS_interaction_method_test := $(filter-out $(BUILD_DIR)/main.o $(BUILD_DIR)/cli/dispatch.o $(BUILD_DIR)/cli/help.o $(BUILD_DIR)/cli/parse.o,$(ALL_OBJS))
TESTOBJS_af_model_test := $(filter-out $(BUILD_DIR)/main.o $(BUILD_DIR)/cli/dispatch.o $(BUILD_DIR)/cli/help.o $(BUILD_DIR)/cli/parse.o,$(ALL_OBJS))
TESTOBJS_partial_spa_test := $(BUILD_DIR)/util/spa_cgf.o
TESTOBJS_method_scope_test := $(filter-out $(BUILD_DIR)/main.o $(BUILD_DIR)/cli/dispatch.o $(BUILD_DIR)/cli/help.o $(BUILD_DIR)/cli/parse.o,$(ALL_OBJS))
TESTOBJS_local_kernel_test := $(filter-out $(BUILD_DIR)/main.o $(BUILD_DIR)/cli/dispatch.o $(BUILD_DIR)/cli/help.o $(BUILD_DIR)/cli/parse.o,$(ALL_OBJS))
TESTOBJS_local_wald_rank_test := $(BUILD_DIR)/interaction/interaction_wald.o $(BUILD_DIR)/util/wald.o $(BUILD_DIR)/util/probability_math.o
TESTOBJS_text_io_test := $(BUILD_DIR)/util/text_stream.o $(ZSTD_OBJS) $(ZLIB_OBJS)
.SECONDEXPANSION:
$(TEST_BUILD_DIR)/%$(EXE): $(TEST_DIR)/%.cpp $$(TESTOBJS_$$*) | tmp
	@mkdir -p $(@D)
	$(CXX) $(TEST_CXXFLAGS) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -MMD -MP -I$(TEST_DIR) $(INCLUDES) $< $(TESTOBJS_$*) -o $@ $(LINK_FLAGS) $(LDLIBS)
-include $(TEST_BINS:$(EXE)=.d)
.PHONY: test run clean
test: $(TEST_BINS)
	@set -e; for t in $(TEST_BINS); do echo "Running $$t"; "$$t"; done
	@echo "ALL TESTS PASSED"
run: all
	$(BIN)
clean:
	rm -rf $(BUILD_DIR) $(DIST_DIR)
# Optional reproducible format-fixture utility (not part of the analysis CLI).
BCF_TOOL := $(BUILD_DIR)/tools/make_bcf_fixture$(EXE)
$(BCF_TOOL): tests/make_bcf_fixture.cpp $(HTSLIB_OBJS) $(HTSCODEC_OBJS) $(ZLIB_OBJS) $(ZSTD_OBJS) $(DEFLATE_OBJS)
	@mkdir -p $(@D)
	$(CXX) $(SPAGXE_CXXFLAGS) $(CPPFLAGS) $(CXXFLAGS) $(INCLUDES) $< $(HTSLIB_OBJS) $(HTSCODEC_OBJS) $(ZLIB_OBJS) $(ZSTD_OBJS) $(DEFLATE_OBJS) -o $@ $(LINK_FLAGS) $(LDLIBS)
.PHONY: fixtures
fixtures: $(BCF_TOOL)
	python3 tests/make_format_fixtures.py
	$(BCF_TOOL) examples/synthetic/synthetic.vcf examples/synthetic/synthetic.bcf

