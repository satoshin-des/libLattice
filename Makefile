CC	= gcc

all: TARGET1 TARGET2

TARGET1: clat.c
	${CC} -Warray-bounds -shared -o libclat.so -c clat.c -lm

TARGET2: main.c
	${CC} -Warray-bounds -o main.exe main.c -L. -lclat -lm
