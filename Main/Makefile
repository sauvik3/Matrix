# This Makefile will build an application using the DLL.

TOOLCHAIN :=	C:\MinGW
SRCTREE :=	C:\Users\admin\Desktop\matrix_cpp\Matrix\Main

CC :=		$(TOOLCHAIN)\bin\g++.exe
WINDRES :=	$(TOOLCHAIN)\bin\windres.exe

SRCDIR :=	$(SRCTREE)\SRC
OUTDIR :=	$(SRCTREE)\BIN
OBJDIR :=	$(SRCTREE)\OBJ
RESDIR :=	$(SRCTREE)\RC
LIBDIR :=	$(SRCTREE)\LIB
INCDIR :=	$(SRCTREE)\INCLUDE

CFLAGS :=	-O3 -std=c++11 -Wall
CFLAGS_EXTRA :=

all: use-static

use-shared: CFLAGS_EXTRA := -DMATRIX_USE_DLL -DDLL_IMPORTS
use-shared: main.exe

use-static: CFLAGS_EXTRA :=
use-static: main.exe

%.o: $(SRCDIR)\%.cpp
	$(CC) -c $(CFLAGS) -I $(INCDIR) $(CFLAGS_EXTRA) -o $(OBJDIR)\$@ $<

resource.o: $(RESDIR)\main.rc
	$(WINDRES) -i $< -o $(OBJDIR)\$@

main.exe: matrix_io.o main.o resource.o
	$(CC) -o $(OUTDIR)\$@ $(patsubst %,$(OBJDIR)\\%,$^) $(LIBDIR)\matrix.lib -s -Wl,-subsystem,console

clean:
	-rm $(OBJDIR)\main.o
	-rm $(OBJDIR)\matrix_io.o
	-rm $(OBJDIR)\resource.o
	-rm $(OUTDIR)\main.exe