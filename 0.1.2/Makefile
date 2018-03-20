OBJS = equation.o method.o runge-kutta.o multi-steps.o ballistic.o
LIBS = -lm `pkg-config --libs gsl`
CFLAGS = -c -O3 -march=native -Wall -Wextra -Wpedantic `pkg-config --cflags gsl`
CC = gcc -g -flto

ballistic: $(OBJS)
	$(CC) $(OBJS) $(LIBS) -o ballistic

equation.o: equation.c equation.h Makefile
	$(CC) $(CFLAGS) equation.c -o equation.o

method.o: method.c method.h equation.h Makefile
	$(CC) $(CFLAGS) method.c -o method.o

runge-kutta.o: runge-kutta.c runge-kutta.h method.h equation.h Makefile
	$(CC) $(CFLAGS) runge-kutta.c -o runge-kutta.o

multi-steps.o: multi-steps.c multi-steps.h runge-kutta.h method.h equation.h \
	Makefile
	$(CC) $(CFLAGS) multi-steps.c -o multi-steps.o

ballistic.o: ballistic.c multi-steps.h runge-kutta.h method.h equation.h \
	Makefile
	$(CC) $(CFLAGS) ballistic.c -o ballistic.o

clean:
	rm -rf html latex *.o ballistic out* *.eps
