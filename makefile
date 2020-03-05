opt = -ffast-math -O3
#cheat = -march=native
ppp = -Wall
a.out: main.o matrixes.o io.o solve.o time.o
	g++ -std=c++17 ${opt} ${ppp} ${cheat} main.o matrixes.o io.o solve.o time.o -pthread  -o  a.out
main.o:  matrixes.h time.h
	g++ -std=c++17 ${opt} ${ppp} ${cheat} -c  main.cpp
matrixes.o:  matrixes.h
	g++ -std=c++17 ${opt} ${ppp} ${cheat} -c matrixes.cpp
io.o: matrixes.h
	g++ -std=c++17 ${opt} ${ppp} ${cheat} -c io.cpp
solve.o: matrixes.h
	g++ -std=c++17 ${opt} ${ppp} ${cheat} -c  solve.cpp
time.o: time.h
	g++ -std=c++17 ${opt} ${ppp} ${cheat} -c  time.cpp
clean:
	rm -f *.o*
