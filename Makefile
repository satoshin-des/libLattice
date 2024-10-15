CC	= gcc

libLattice: clat.c
	${CC} -shared -o libclat.so -c clat.c -lm

main: main.c
	${CC} -o main.exe main.c -L. -lclat -lm

NTL:
	sudo apt-get install -y libntl-dev