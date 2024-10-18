CC	= gcc

libLattice: clat.c
	${CC} -shared -o libclat.so -c clat.c -lm

main: main.c
	${CC} -o main.exe main.c -L. -lclat -lm

# If you have not installed NTL library(Victor Shoup, https://libntl.org/), please execute the next command and install it.
NTL:
	sudo apt-get install -y libntl-dev