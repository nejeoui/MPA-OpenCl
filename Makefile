# MPA-OpenCL
#
#   make            build everything
#   make test       build and run the GMP-vs-device correctness harness
#   make bench      build the three timing hosts
#
# Dependencies: an OpenCL 1.2+ ICD, GMP, and (for the benchmark hosts only)
# OpenSSL.
#
#   Debian/Ubuntu:  sudo apt install ocl-icd-opencl-dev libgmp-dev libssl-dev
#                   plus a runtime, e.g. pocl-opencl-icd, or your vendor's driver
#   macOS:          brew install gmp openssl
#                   (OpenCL ships with the OS; Homebrew paths are found below)

CC      ?= cc
CFLAGS  ?= -O2 -std=c99 -Wall -Wextra -Wno-unused-parameter
CPPFLAGS += -DCL_TARGET_OPENCL_VERSION=120 -D_POSIX_C_SOURCE=200809L

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
  CL_LDLIBS := -framework OpenCL
  # _POSIX_C_SOURCE hides Darwin extensions such as _SC_NPROCESSORS_ONLN;
  # _DARWIN_C_SOURCE puts them back.
  CPPFLAGS += -D_DARWIN_C_SOURCE
  # Homebrew does not install into the default search path, and openssl is
  # keg-only, so ask brew where each package actually lives. Works on both
  # /opt/homebrew (Apple Silicon) and /usr/local (Intel). If brew is absent
  # these expand to nothing and the build falls back to the system paths.
  GMP_PREFIX := $(shell brew --prefix gmp 2>/dev/null)
  SSL_PREFIX := $(shell brew --prefix openssl@3 2>/dev/null || brew --prefix openssl 2>/dev/null)
  ifneq ($(GMP_PREFIX),)
    CPPFLAGS += -I$(GMP_PREFIX)/include
    LDFLAGS  += -L$(GMP_PREFIX)/lib
  endif
  ifneq ($(SSL_PREFIX),)
    CPPFLAGS += -I$(SSL_PREFIX)/include
    LDFLAGS  += -L$(SSL_PREFIX)/lib
  endif
else
  CL_LDLIBS := -lOpenCL
endif

BENCH := mpa_8bits mpa_16bits mpa_32bits

.PHONY: all test bench perf compare ecdsa run clean

all: mpa_test mpa_bench mpa_compare ecdsa_bench mpa_run $(BENCH)

# The dependency-free host: C99 and an OpenCL ICD, nothing else. The explicit
# rule keeps it away from the mpa_% pattern below, which links GMP and OpenSSL.
mpa_run: mpa_run.c
	$(CC) $(CFLAGS) $(CPPFLAGS) $< -o $@ $(LDFLAGS) $(CL_LDLIBS)

mpa_test: mpa_test.c mpa_ref.h
	$(CC) $(CFLAGS) $(CPPFLAGS) $< -o $@ $(LDFLAGS) $(CL_LDLIBS) -lgmp

mpa_bench: mpa_bench.c mpa_ref.h
	$(CC) $(CFLAGS) $(CPPFLAGS) $< -o $@ $(LDFLAGS) $(CL_LDLIBS) -lgmp

# OpenMP. Without it the multi-threaded rows collapse to one thread and the
# CPU baseline is understated by the core count, so mpa_compare warns loudly
# at run time if that happens. Apple clang ships no OpenMP runtime: it needs
# `brew install libomp` and the -Xpreprocessor spelling below.
ifeq ($(UNAME_S),Darwin)
  OMP_PREFIX := $(shell brew --prefix libomp 2>/dev/null)
  ifneq ($(OMP_PREFIX),)
    OMPFLAGS := -Xpreprocessor -fopenmp -I$(OMP_PREFIX)/include
    OMPLIBS  := -L$(OMP_PREFIX)/lib -lomp
  endif
else
  OMPFLAGS := $(shell $(CC) -fopenmp -E -x c /dev/null >/dev/null 2>&1 && echo -fopenmp)
  OMPLIBS  :=
endif

ecdsa_bench: ecdsa_bench.c ecdsaKernel.cl
	$(CC) $(CFLAGS) $(OMPFLAGS) $(CPPFLAGS) -Wno-deprecated-declarations $< -o $@ $(LDFLAGS) $(OMPLIBS) $(CL_LDLIBS) -lgmp -lcrypto

mpa_compare: mpa_compare.c mpa_ref.h
	$(CC) $(CFLAGS) $(OMPFLAGS) $(CPPFLAGS) $< -o $@ $(LDFLAGS) $(OMPLIBS) $(CL_LDLIBS) -lgmp -lcrypto

# The kernels are read from the working directory at run time, and their
# #include of the matching .h is resolved with -I. inside the harness.
test: mpa_test
	./mpa_test $(TESTARGS)

bench: $(BENCH)

# Performance sweep: every optimization level, verified before it is timed.
perf: mpa_bench
	./mpa_bench $(PERFARGS)

# Head-to-head against GMP and OpenSSL on identical, verified operands.
compare: mpa_compare
	./mpa_compare $(CMPARGS)

# Minimal host: drive the kernels with no GMP and no OpenSSL.
run: mpa_run
	./mpa_run $(RUNARGS)

# Batched P-256 ECDSA verification, validated against OpenSSL.
ecdsa: ecdsa_bench
	./ecdsa_bench $(ECDSAARGS)

mpa_%: mpa_%.c
	$(CC) $(CFLAGS) $(CPPFLAGS) $< -o $@ $(LDFLAGS) $(CL_LDLIBS) -lgmp -lcrypto

clean:
	rm -f mpa_test mpa_bench mpa_compare ecdsa_bench mpa_run $(BENCH)
