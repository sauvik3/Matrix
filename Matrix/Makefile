# This Makefile will build the Matrix Library.

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


dirs:
	@if [ ! -d $(OUTDIR) ]; then mkdir -p $(notdir $(OUTDIR)); fi
	@if [ ! -d $(OBJDIR) ]; then mkdir -p $(notdir $(OBJDIR)); fi
	@if [ ! -d $(LIBDIR) ]; then mkdir -p $(notdir $(LIBDIR)); fi

all: static

dll: CFLAGS_EXTRA := -fpic -DMATRIX_USE_DLL -DDLL_EXPORTS
dll: dirs Matrix.so

static: CFLAGS_EXTRA :=
static: dirs Matrix.a

Matrix.so: matrix_exception.o matrix.o

Matrix.a: matrix_exception.o matrix.o

%.o: $(SRCDIR)/%.cpp
	$(CXX) -c $(CFLAGS) -I $(INCDIR) $(CFLAGS_EXTRA) -o $(OBJDIR)/$@ $<

%.so: 
	$(CXX) -shared -s -o $(OUTDIR)/$@ $(patsubst %,$(OBJDIR)/%,$^)

%.a: 
	$(AR) -ru $(LIBDIR)/$@ $(patsubst %,$(OBJDIR)/%,$^)


clean::
	-rm $(OBJDIR)/matrix.o
	-rm $(LIBDIR)/Matrix.a
	-rm $(OUTDIR)/Matrix.so