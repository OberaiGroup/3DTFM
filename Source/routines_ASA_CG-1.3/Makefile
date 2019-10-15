default: all

CC = gcc -lm -W -Wall -pedantic -Wmissing-prototypes \
	-Wredundant-decls -Wnested-externs -Wdisabled-optimization \
	-ansi -g -fexceptions -Wno-parentheses -Wshadow -Wcast-align \
	-Winline -Wstrict-prototypes -Wno-unknown-pragmas

# for speed:
#CC = gcc -lm -O3

all: driverMain

driverMain: Makefile driverMain.c asa_cg.o driverASA_CG.o driverMain2.o
	$(CC) driverMain.c -o driverMain asa_cg.o driverASA_CG.o driverMain2.o

asa_cg.o: asa_cg.c
	$(CC) -c $*.c
driverASA_CG.o: driverASA_CG.c
	$(CC) -c $*.c
driverMain2.o: driverMain2.c
	$(CC) -c $*.c

clean:
	rm *.o

purge:
	rm *.o driverMain
