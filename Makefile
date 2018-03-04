# This Makefile will build the Project recursively.


all: Bench-static Main-static

Matrix-static: 
		$(MAKE) -C Matrix static

Bench-static:	Matrix-static
		mkdir -p Bench/lib
		cp -p Matrix/lib/Matrix.a -t Bench/lib
		$(MAKE) -C Bench use-static

Main-static:	Matrix-static
		mkdir -p Main/lib
		cp -p Matrix/lib/Matrix.a -t Main/lib
		$(MAKE) -C Main use-static
		
Matrix-dll: 
		$(MAKE) -C Matrix dll

Bench-dll:	Matrix-dll
		mkdir -p Bench/bin
		cp -p Matrix/bin/Matrix.so -t Bench/bin
		$(MAKE) -C Bench use-shared

Main-dll:	Matrix-dll
		mkdir -p Main/bin
		cp -p Matrix/bin/Matrix.so -t Main/bin
		$(MAKE) -C Main use-shared

clean: 
		$(MAKE) -C Bench	clean
		$(MAKE) -C Main	clean
		$(MAKE) -C Matrix	clean
		
.PHONY: all Matrix clean
