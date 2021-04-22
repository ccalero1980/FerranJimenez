.SUFFIXES:
.SUFFIXES: .c .o
.c.o:
	gcc -c Barra_fluid.c

ex_barra_fluid: Barra_fluid.o
	gcc -o ex_barra_fluid Barra_fluid.o -lm

clean:
	rm *.o
