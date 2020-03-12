opt = -ffast-math -O3
#cheat = -march=native
ppp = -Wall
a.out: main.o matrixes.o io.o solve.o
	mpicxx -std=c++17 ${opt} ${ppp} ${cheat} main.o matrixes.o io.o solve.o -o a.out
main.o:  matrixes.h
	mpicxx -std=c++17 ${opt} ${ppp} ${cheat} -c  main.cpp
matrixes.o:  matrixes.h
	mpicxx -std=c++17 ${opt} ${ppp} ${cheat} -c matrixes.cpp
io.o: matrixes.h
	mpicxx -std=c++17 ${opt} ${ppp} ${cheat} -c io.cpp
solve.o: matrixes.h
	mpicxx -std=c++17 ${opt} ${ppp} ${cheat} -c  solve.cpp
clean:
	rm -f *.o*
