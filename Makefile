CFLAGS = -O3 -std=c99 
LDFLAGS = -lm
CC = gcc

relax : relax.c
	gcc -o relax relax.c -lm $(CFLAGS)
