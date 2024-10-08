CC	= gcc

all: TARGET1 TARGET2

TARGET1: clat.c
	${CC} -shared -o libclat.so -c clat.c -lm

TARGET2: main.c
	${CC} -o main main.c -L. -lclat -lm
