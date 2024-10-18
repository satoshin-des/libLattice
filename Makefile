CC	= gcc
SRC = hoge.c

libLattice: clat.c
	${CC} -shared -o libclat.so -c clat.c -lm

main: main.c
	${CC} -o main.exe main.c -L. -lclat -lm

hoge: ${SRC}
	${CC} ${SRC} -L. -lclat -lm
