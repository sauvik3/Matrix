# This Makefile will build an application using the DLL.

SRCTREE :=	./

CXX ?=		g++
AR ?=		ar

SRCDIR :=	$(SRCTREE)src
OUTDIR :=	$(SRCTREE)bin
OBJDIR :=	$(SRCTREE)obj
LIBDIR :=	$(SRCTREE)lib
INCDIR :=	$(SRCTREE)include

CFLAGS :=	-O3 -std=c++11 -Wall -pedantic
CFLAGS_EXTRA :=
LFLAGS :=


dirs:
	@if [ ! -d $(OUTDIR) ]; then mkdir -p $(notdir $(OUTDIR)); fi
	@if [ ! -d $(OBJDIR) ]; then mkdir -p $(notdir $(OBJDIR)); fi
	@if [ ! -d $(LIBDIR) ]; then mkdir -p $(notdir $(LIBDIR)); fi

all: use-static

use-shared: CFLAGS_EXTRA := -pthread -DMATRIX_USE_DLL -DDLL_IMPORTS
use-shared: LFLAGS := -L$(SRCTREE) -l:$(OUTDIR)/Matrix.so
use-shared: dirs main

use-static: CFLAGS_EXTRA := -pthread
use-static: LFLAGS := $(LIBDIR)/Matrix.a
use-static: dirs main

%.o: $(SRCDIR)/%.cpp
	$(CXX) -c $(CFLAGS) -I $(INCDIR) $(CFLAGS_EXTRA) -o $(OBJDIR)/$@ $<

main: matrix_io.o main.o
	$(CXX) -o $(OUTDIR)/$@ $(patsubst %,$(OBJDIR)/%,$^) $(LFLAGS) -s -pthread

clean:
	-rm $(OBJDIR)/main.o
	-rm $(OBJDIR)/matrix_io.o
	-rm $(OUTDIR)/main
	-rm $(LIBDIR)/Matrix.a
	-rm $(OUTDIR)/Matrix.so