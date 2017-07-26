# This Makefile will build the Matrix Library.

SRCTREE :=	.

CXX ?=		g++
AR ?=		ar

SRCDIR :=	$(SRCTREE)/src
OUTDIR :=	$(SRCTREE)/bin
OBJDIR :=	$(SRCTREE)/obj
LIBDIR :=	$(SRCTREE)/lib
INCDIR :=	$(SRCTREE)/include

CFLAGS :=	-O3 -std=c++11 -fpic -Wall -pedantic
CFLAGS_EXTRA :=


dirs:
	@if [ ! -d $(OUTDIR) ]; then mkdir -p $(notdir $(OUTDIR)); fi
	@if [ ! -d $(OBJDIR) ]; then mkdir -p $(notdir $(OBJDIR)); fi
	@if [ ! -d $(LIBDIR) ]; then mkdir -p $(notdir $(LIBDIR)); fi

all: static

dll: CFLAGS_EXTRA := -DMATRIX_USE_DLL -DDLL_EXPORTS
dll: dirs matrix.so

static: CFLAGS_EXTRA :=
static: dirs matrix.a

matrix.so: matrix_exception.o matrix.o

matrix.a: matrix_exception.o matrix.o

%.o: $(SRCDIR)/%.cpp
	$(CXX) -c $(CFLAGS) -I $(INCDIR) $(CFLAGS_EXTRA) -o $(OBJDIR)/$@ $<

%.so: 
	$(CXX) -shared -s -o $(OUTDIR)/$@ $(patsubst %,$(OBJDIR)/%,$^)

%.a: 
	$(AR) -ru $(LIBDIR)/$@ $(patsubst %,$(OBJDIR)/%,$^)


clean::
	-rm $(OBJDIR)/matrix.o
	-rm $(LIBDIR)/matrix.a
	-rm $(OUTDIR)/matrix.so
