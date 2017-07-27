# This Makefile will build the Project recursively.


all: Bench-static Main-static

Matrix-static: 
		make -C Matrix static

Bench-static:	Matrix-static
		mkdir -p Bench/lib
		cp -p Matrix/lib/Matrix.a -t Bench/lib
		make -C Bench use-static

Main-static:	Matrix-static
		mkdir -p Main/lib
		cp -p Matrix/lib/Matrix.a -t Main/lib
		make -C Main use-static
		
Matrix-dll: 
		make -C Matrix dll

Bench-dll:	Matrix-dll
		mkdir -p Bench/bin
		cp -p Matrix/bin/libMatrix.so -t Bench/bin
		make -C Bench use-shared

Main-dll:	Matrix-dll
		mkdir -p Main/bin
		cp -p Matrix/bin/libMatrix.so -t Main/bin
		make -C Main use-shared

clean: 
		make -C Bench	clean
		make -C Main	clean
		make -C Matrix	clean
		
.PHONY: all Matrix clean
