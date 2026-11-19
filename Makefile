# SPDX-FileCopyrightText: 2020-2026 The psvpp-ba authors (see AUTHORS)
# SPDX-License-Identifier: GPL-3.0-or-later
# psvpp-ba - build system
#
#   make                 optimised build (Gurobi is used when found, see below)
#   make debug           unoptimised build with debug symbols and sanitizers
#   make test            build and run the test suite
#   make install         copy the program to $(PREFIX)/bin (PREFIX=/usr/local)
#   make clean           remove build artefacts
#   make help            show this help
#
# Gurobi: the branch-and-cut (`bc`) and complete-model (`complete`) commands
# need Gurobi. It is taken from GUROBI_HOME (set by the Gurobi installer or an
# environment module), e.g. GUROBI_HOME=/opt/gurobi1400/linux64; when that is
# unset, or names a folder without include/gurobi_c.h (a removed version), the
# newest installation in /opt or /Library is used. Use `make GUROBI=0` to build
# without it (the heuristic and `routes` still work).
#
# OpenMP: used to run independent heuristic runs in parallel (`--threads`).
# It is used when the compiler supports it; disable with `make OPENMP=0`.

# ---------------------------------------------------------------- toolchain
CC      ?= cc
BUILD   ?= build
PREFIX  ?= /usr/local
TARGET  := $(BUILD)/psvpp-ba
VERSION := 1.0.0

# Newest C standard supported by the compiler (C23, or C2x on GCC 13 / Clang 16-17).
STD := $(shell $(CC) -std=c23 -E -x c /dev/null >/dev/null 2>&1 && echo c23 || echo c2x)

WARNINGS := -Wall -Wextra -Wpedantic -Wshadow -Wformat=2 -Wstrict-prototypes \
            -Wmissing-prototypes -Wold-style-definition -Wvla -Wcast-qual -Wundef
# CFLAGS/LDFLAGS are for the user (e.g. make CFLAGS="-O3 -march=native");
# the flags the build needs are kept in separate variables.
CFLAGS   ?= -O2
LDFLAGS  ?=
# Strict floating-point arithmetic, so that the costs do not depend on the compiler or the CPU:
# -ffp-contract=off forbids fused multiply-add, and -fno-fast-math switches off the relaxed
# arithmetic that some compilers enable by default (Intel icx). The pruning of insertion
# evaluations compares sums of floats exactly and relies on it (tests/test_moves.c checks this).
ALL_CFLAGS   = -std=$(STD) $(WARNINGS) -ffp-contract=off -fno-fast-math -MMD -MP $(OMP_FLAGS) $(CFLAGS)
ALL_CPPFLAGS = -Isrc -DPSVPP_VERSION='"$(VERSION)"' $(GRB_CPPFLAGS) $(CPPFLAGS)
ALL_LDFLAGS  = $(OMP_FLAGS) $(GRB_LDFLAGS) $(LDFLAGS)
ALL_LDLIBS   = $(GRB_LDLIBS) -lm $(LDLIBS)

# ---------------------------------------------------------------- OpenMP
OPENMP ?= auto
ifeq ($(OPENMP),auto)
  OPENMP := $(shell echo 'int main(void){return 0;}' | $(CC) -fopenmp -x c - -o /dev/null >/dev/null 2>&1 && echo 1 || echo 0)
endif
ifeq ($(OPENMP),1)
  OMP_FLAGS := -fopenmp
endif

# ---------------------------------------------------------------- Gurobi
GUROBI ?= auto
ifneq ($(GUROBI),0)
  # newest installation in the standard places, by version number (gurobi1400 is newer than gurobi951)
  GUROBI_FOUND := $(shell for d in /opt/gurobi*/linux64 /Library/gurobi*/macos_universal2; do \
                    [ -f "$$d/include/gurobi_c.h" ] && echo "$$(echo "$$d" | sed 's|.*/gurobi\([0-9]*\)/.*|\1|') $$d"; \
                  done | sort -n | tail -n 1 | cut -d' ' -f2)
  ifeq ($(GUROBI_HOME),)
    GUROBI_HOME := $(GUROBI_FOUND)
  else ifeq ($(wildcard $(GUROBI_HOME)/include/gurobi_c.h),)
    $(warning GUROBI_HOME=$(GUROBI_HOME) has no include/gurobi_c.h; using $(or $(GUROBI_FOUND),no Gurobi))
    GUROBI_HOME := $(GUROBI_FOUND)
  endif
  GUROBI_LIB := $(firstword $(patsubst lib%.so,%,$(notdir $(wildcard \
                  $(GUROBI_HOME)/lib/libgurobi[0-9][0-9].so $(GUROBI_HOME)/lib/libgurobi[0-9][0-9][0-9].so \
                  $(GUROBI_HOME)/lib/libgurobi[0-9][0-9].dylib $(GUROBI_HOME)/lib/libgurobi[0-9][0-9][0-9].dylib))))
  GUROBI_LIB := $(patsubst lib%.dylib,%,$(GUROBI_LIB))
  ifneq ($(GUROBI_LIB),)
    GRB_CPPFLAGS := -DPSVPP_WITH_GUROBI -isystem $(GUROBI_HOME)/include
    GRB_LDFLAGS  := -L$(GUROBI_HOME)/lib -Wl,-rpath,$(GUROBI_HOME)/lib
    GRB_LDLIBS   := -l$(GUROBI_LIB)
    WITH_GUROBI := yes
  else ifeq ($(GUROBI),1)
    $(error Gurobi requested but not found; set GUROBI_HOME, e.g. GUROBI_HOME=/opt/gurobi1400/linux64)
  endif
endif

# ---------------------------------------------------------------- sources
COMMON_SRC    := $(wildcard src/common/*.c)
HEURISTIC_SRC := $(wildcard src/heuristic/*.c)
# the two MIP models need Gurobi; the voyage enumeration and the scheduling algorithms do not
GUROBI_SRC    := src/exact/grb.c src/exact/reduced_model.c src/exact/complete_model.c \
                 src/exact/subproblem_models.c
EXACT_SRC     := $(filter-out $(GUROBI_SRC),$(wildcard src/exact/*.c))
ifeq ($(WITH_GUROBI),yes)
  EXACT_SRC += $(GUROBI_SRC)
endif
LIB_SRC := $(COMMON_SRC) $(HEURISTIC_SRC) $(EXACT_SRC)
SRC     := src/main.c $(LIB_SRC)
OBJ     := $(SRC:%.c=$(BUILD)/%.o)
LIB_OBJ := $(LIB_SRC:%.c=$(BUILD)/%.o)

TEST_SRC := $(wildcard tests/*.c)
TEST_BIN := $(TEST_SRC:tests/%.c=$(BUILD)/tests/%)

# ---------------------------------------------------------------- rules
.PHONY: all debug test install clean help format FORCE

all: $(TARGET)
	@echo "built $(TARGET) (C standard: $(STD), Gurobi: $(if $(WITH_GUROBI),$(GUROBI_LIB),no), OpenMP: $(if $(filter 1,$(OPENMP)),yes,no))"

$(TARGET): $(OBJ)
	$(CC) $(ALL_LDFLAGS) -o $@ $^ $(ALL_LDLIBS)

# Objects are rebuilt when the compilation flags change (e.g. GUROBI=0 vs GUROBI=1).
FLAGS_STAMP := $(BUILD)/.flags
$(FLAGS_STAMP): FORCE
	@mkdir -p $(BUILD)
	@echo '$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) $(ALL_LDFLAGS) $(ALL_LDLIBS)' | cmp -s - $@ || \
	  echo '$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) $(ALL_LDFLAGS) $(ALL_LDLIBS)' > $@

$(BUILD)/%.o: %.c $(FLAGS_STAMP)
	@mkdir -p $(dir $@)
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) -c -o $@ $<

debug:
	$(MAKE) BUILD=build-debug CFLAGS="-O0 -g3 -fsanitize=address,undefined -fno-omit-frame-pointer" \
	        LDFLAGS="-fsanitize=address,undefined"

$(BUILD)/tests/%: tests/%.c $(LIB_OBJ) $(FLAGS_STAMP)
	@mkdir -p $(dir $@)
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) $(ALL_LDFLAGS) -o $@ $< $(LIB_OBJ) $(ALL_LDLIBS)

test: $(TARGET) $(TEST_BIN)
	@for t in $(TEST_BIN); do echo "== $$t"; $$t || exit 1; done
	@./tests/smoke.sh $(TARGET)

install: $(TARGET)
	install -d $(DESTDIR)$(PREFIX)/bin
	install -m 755 $(TARGET) $(DESTDIR)$(PREFIX)/bin/psvpp-ba

format:
	clang-format -i src/*.c src/*/*.c src/*/*.h tests/*.c tests/*.h

clean:
	rm -rf build build-debug

help:
	@sed -n '2,/^# ---/p' Makefile | sed '$$d' | sed 's/^# \{0,1\}//'

-include $(OBJ:.o=.d)
