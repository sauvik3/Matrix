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
CFLAGS_EXTRA :=  -coverage
LFLAGS :=	-lgcov --coverage


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

use-shared: CFLAGS_EXTRA += -pthread -DMATRIX_USE_DLL -DDLL_IMPORTS
use-shared: LFLAGS += -L$(SRCTREE) -l:$(OUTDIR)/Matrix.so
use-shared: dirs main

use-static: CFLAGS_EXTRA += -pthread
use-static: LFLAGS += $(LIBDIR)/Matrix.a
use-static: dirs main

%.o: $(SRCDIR)/%.cpp
	$(CXX) -c $(CFLAGS) -I $(INCDIR) $(CFLAGS_EXTRA) -o $(OBJDIR)/$@ $<

main: main.o
	$(CXX) -o $(OUTDIR)/$@ $(patsubst %,$(OBJDIR)/%,$^) $(LFLAGS) -s -pthread

.PHONY : clean
clean:
	$(RM) $(OBJDIR)/main.o
	$(RM) $(OUTDIR)/main
	$(RM) $(LIBDIR)/Matrix.a
	$(RM) $(OUTDIR)/Matrix.so	
