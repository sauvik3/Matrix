# This Makefile will build the Matrix Library.

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

all: static

dll: CFLAGS_EXTRA := -pthread -fpic -DMATRIX_USE_DLL -DDLL_EXPORTS
dll: dirs Matrix.so

static: CFLAGS_EXTRA := -pthread
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
	$(RM) $(OBJDIR)/matrix.o
	$(RM) $(LIBDIR)/Matrix.a
	$(RM) $(OUTDIR)/Matrix.so
