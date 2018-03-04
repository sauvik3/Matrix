# This Makefile will build an application using the DLL.

SRCTREE :=	./

CXX ?=		g++
AR ?=		ar

SRCDIR :=	$(SRCTREE)src
OUTDIR :=	$(SRCTREE)bin
OBJDIR :=	$(SRCTREE)obj
LIBDIR :=	$(SRCTREE)lib
INCDIR :=	$(SRCTREE)include

CFLAGS :=	-O3 -std=c++11 -Wall -pedantic -coverage
CFLAGS_EXTRA :=
LFLAGS :=


dirs:
ifeq ("$(wildcard $(OBJDIR)/.)", "")
	-mkdir $(notdir $(OBJDIR))
endif
ifeq ("$(wildcard $(LIBDIR)/.)", "")
	-mkdir $(notdir $(LIBDIR))
endif
ifeq ("$(wildcard $(OUTDIR)/.)", "")
	-mkdir $(notdir $(OUTDIR))
endif

all: use-static

use-shared: CFLAGS_EXTRA := -pthread -DMATRIX_USE_DLL -DDLL_IMPORTS
use-shared: LFLAGS := -L$(SRCTREE) -l:$(OUTDIR)/Matrix.so -lgcov --coverage
use-shared: dirs main

use-static: CFLAGS_EXTRA := -pthread
use-static: LFLAGS := $(LIBDIR)/Matrix.a -lgcov --coverage
use-static: dirs main

%.o: $(SRCDIR)/%.cpp
	$(CXX) -c $(CFLAGS) -I $(INCDIR) $(CFLAGS_EXTRA) -o $(OBJDIR)/$@ $<

main: matrix_io.o main.o
	$(CXX) -o $(OUTDIR)/$@ $(patsubst %,$(OBJDIR)/%,$^) $(LFLAGS) -s -pthread

clean:
	$(RM) $(OBJDIR)/main.o
	$(RM) $(OBJDIR)/matrix_io.o
	$(RM) $(OUTDIR)/main
	$(RM) $(LIBDIR)/Matrix.a
	$(RM) $(OUTDIR)/Matrix.so
