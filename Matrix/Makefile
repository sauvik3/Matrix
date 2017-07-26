# This Makefile will build the Matrix Library.

SRCTREE :=	.

CXX ?=		g++
AR ?=		ar

SRCDIR :=	$(SRCTREE)\src
OUTDIR :=	$(SRCTREE)\bin
OBJDIR :=	$(SRCTREE)\obj
RESDIR :=	$(SRCTREE)\rc
LIBDIR :=	$(SRCTREE)\lib
INCDIR :=	$(SRCTREE)\include

CFLAGS :=	-O3 -std=c++11 -fpic -Wall -pedantic
CFLAGS_EXTRA :=


all: static

dll: CFLAGS_EXTRA := -DMATRIX_USE_DLL -DDLL_EXPORTS
dll: matrix.so

static: CFLAGS_EXTRA :=
static: matrix.a

matrix.so: matrix_exception.o matrix.o

matrix.a: matrix_exception.o matrix.o

%.o: $(SRCDIR)\%.cpp
	$(CXX) -c $(CFLAGS) -I $(INCDIR) $(CFLAGS_EXTRA) -o $(OBJDIR)\$@ $<

%.so: 
	$(CXX) -shared -s -o $(OUTDIR)\$@ $(patsubst %,$(OBJDIR)\\%,$^)

%.a: 
	$(AR) -ru $(OUTDIR)\$@ $(patsubst %,$(OBJDIR)\\%,$^)


clean::
	-rm $(OBJDIR)\matrix.o
	-rm $(OUTDIR)\matrix.a
	-rm $(OUTDIR)\matrix.so
