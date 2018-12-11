.PHONY: clean strip

PGOOBJS = utils.pgo equation.pgo method.pgo runge-kutta.pgo multi-steps.pgo \
	ballistic.pgo
OBJS = utils.o equation.o method.o runge-kutta.o multi-steps.o ballistic.o
LIBS = -lm `pkg-config --libs gsl libxml-2.0 glib-2.0`
CFLAGS = `pkg-config --cflags gsl libxml-2.0 glib-2.0` -c -O3 -march=native \
	-Wall -Wextra -Wpedantic -D_FORTIFY_SOURCE=2
PGOGEN = -fprofile-generate
PGOUSE = -fprofile-use -fprofile-correction
CC = gcc -g -flto

ballistic: $(OBJS)
	$(CC) $(OBJS) $(LIBS) -o ballistic

ballisticpgo: $(PGOOBJS)
	$(CC) $(PGOGEN) $(PGOOBJS) $(LIBS) -o ballisticpgo
	./ballisticpgo ../tests/case0.xml out
	./ballisticpgo ../tests/case1.xml out
	./ballisticpgo ../tests/case2.xml out
	./ballisticpgo ../tests/case2ms.xml out
	./ballisticpgo ../tests/case3.xml out
	./ballisticpgo ../tests/case-s.xml out
	sh scriptpgo.sh

utils.pgo: utils.c utils.h config.h Makefile
	$(CC) $(CFLAGS) $(PGOGEN) utils.c -o utils.pgo

equation.pgo: equation.c equation.h utils.h config.h Makefile
	$(CC) $(CFLAGS) $(PGOGEN) equation.c -o equation.pgo

method.pgo: method.c method.h equation.h utils.h config.h Makefile
	$(CC) $(CFLAGS) $(PGOGEN) method.c -o method.pgo

runge-kutta.pgo: runge-kutta.c runge-kutta.h method.h equation.h utils.h \
	config.h Makefile
	$(CC) $(CFLAGS) $(PGOGEN) runge-kutta.c -o runge-kutta.pgo

multi-steps.pgo: multi-steps.c multi-steps.h runge-kutta.h method.h equation.h \
	utils.h config.h Makefile
	$(CC) $(CFLAGS) $(PGOGEN) multi-steps.c -o multi-steps.pgo

ballistic.pgo: ballistic.c multi-steps.h runge-kutta.h method.h equation.h \
	utils.h config.h Makefile
	$(CC) $(CFLAGS) $(PGOGEN) ballistic.c -o ballistic.pgo

utils.o: ballisticpgo utils.gcda
	$(CC) $(CFLAGS) $(PGOUSE) utils.c -o utils.o

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
