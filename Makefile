CC	= gcc


libLattice: clat.c
	${CC} -shared -o libclat.so -c clat.c -lm

main: main.c
	${CC} -o main.exe main.c -L. -lclat -lm
