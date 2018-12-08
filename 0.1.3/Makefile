.PHONY: clean strip

PGOOBJS = equation.pgo method.pgo runge-kutta.pgo multi-steps.pgo ballistic.pgo
OBJS = equation.o method.o runge-kutta.o multi-steps.o ballistic.o
LIBS = -lm `pkg-config --libs gsl`
CFLAGS = `pkg-config --cflags gsl` -c -O3 -march=native -Wall -Wextra \
	-Wpedantic -D_FORTIFY_SOURCE=2
PGOGEN = -fprofile-generate
PGOUSE = -fprofile-use -fprofile-correction
CC = gcc -g -flto

ballistic: $(OBJS)
	$(CC) $(OBJS) $(LIBS) -o ballistic

ballisticpgo: $(PGOOBJS)
	$(CC) $(PGOGEN) $(PGOOBJS) $(LIBS) -o ballisticpgo
	cp script.sh scriptpgo.sh
	sed -i 's/ballistic/ballisticpgo/g' scriptpgo.sh
	sh scriptpgo.sh

equation.pgo: equation.c equation.h Makefile
	$(CC) $(CFLAGS) $(PGOGEN) equation.c -o equation.pgo

method.pgo: method.c method.h equation.h Makefile
	$(CC) $(CFLAGS) $(PGOGEN) method.c -o method.pgo

runge-kutta.pgo: runge-kutta.c runge-kutta.h method.h equation.h Makefile
	$(CC) $(CFLAGS) $(PGOGEN) runge-kutta.c -o runge-kutta.pgo

multi-steps.pgo: multi-steps.c multi-steps.h runge-kutta.h method.h equation.h \
	Makefile
	$(CC) $(CFLAGS) $(PGOGEN) multi-steps.c -o multi-steps.pgo

ballistic.pgo: ballistic.c multi-steps.h runge-kutta.h method.h equation.h \
	Makefile
	$(CC) $(CFLAGS) $(PGOGEN) ballistic.c -o ballistic.pgo

equation.o: ballisticpgo equation.gcda
	$(CC) $(CFLAGS) $(PGOUSE) equation.c -o equation.o

method.o: ballisticpgo method.gcda
	$(CC) $(CFLAGS) $(PGOUSE) method.c -o method.o

runge-kutta.o: ballisticpgo runge-kutta.gcda 
	$(CC) $(CFLAGS) $(PGOUSE) runge-kutta.c -o runge-kutta.o

multi-steps.o: ballisticpgo multi-steps.gcda 
	$(CC) $(CFLAGS) $(PGOUSE) multi-steps.c -o multi-steps.o

ballistic.o: ballisticpgo ballistic.gcda 
	$(CC) $(CFLAGS) $(PGOUSE) ballistic.c -o ballistic.o

clean:
	rm -rf html latex *pgo *.o ballistic out* *.eps

strip:
	strip ballistic
