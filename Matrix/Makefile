# This Makefile will build the Matrix Library.

TOOLCHAIN :=	C:\MinGW
SRCTREE :=	C:\Users\admin\Desktop\matrix_cpp\Matrix\Matrix

CC :=		$(TOOLCHAIN)\bin\g++.exe
AR :=		$(TOOLCHAIN)\bin\ar.exe
WINDRES :=	$(TOOLCHAIN)\bin\windres.exe

SRCDIR :=	$(SRCTREE)\SRC
OUTDIR :=	$(SRCTREE)\BIN
OBJDIR :=	$(SRCTREE)\OBJ
RESDIR :=	$(SRCTREE)\RC
LIBDIR :=	$(SRCTREE)\LIB
INCDIR :=	$(SRCTREE)\INCLUDE

CFLAGS :=	-O3 -std=c++11 -Wall
CFLAGS_EXTRA :=


all: static

dll: CFLAGS_EXTRA := -DMATRIX_USE_DLL -DDLL_EXPORTS
dll: matrix.dll

static: CFLAGS_EXTRA :=
static: matrix.lib

matrix.dll: matrix_exception.o matrix.o resource.o

matrix.lib: matrix_exception.o matrix.o

%.o: $(SRCDIR)\%.cpp
	$(CC) -c $(CFLAGS) -I $(INCDIR) $(CFLAGS_EXTRA) -o $(OBJDIR)\$@ $<

%.dll: 
	$(CC) -shared -s -o $(OUTDIR)\$@ $(patsubst %,$(OBJDIR)\\%,$^) -Wl,--subsystem,windows,--output-def,$(OUTDIR)\$*.def,--out-implib,$(OUTDIR)\$*.lib

%.lib: 
	$(AR) -ru $(OUTDIR)\$@ $(patsubst %,$(OBJDIR)\\%,$^)

resource.o: $(RESDIR)\matrix.rc
	$(WINDRES) -i $< -o $(OBJDIR)\$@

clean::
	-rm $(OBJDIR)\matrix.o
	-rm $(OBJDIR)\resource.o
	-rm $(OUTDIR)\matrix.def
	-rm $(OUTDIR)\matrix.lib
	-rm $(OUTDIR)\matrix.dll