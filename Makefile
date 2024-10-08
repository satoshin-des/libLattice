CC	= gcc
SRC	= clat.c
CFLAGS	= -shared -o libclat.so -c
LDFLAGS = -lm

all:
	${CC} ${CFLAGS} ${SRC} ${LDFLAGS}
