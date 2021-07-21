# EXTRA lets you do:  make -DEXTRA=blahblah
# If you want to compile for debugging, change "-O" to "-g"
# To turn on profiling, uncomment next line

CC=gcc

#PROFILE = -pg
# gsl info page recommends -ansi and -pedantic.  -ansi gives errors,
# and -pedantic gives warnings (e.g. about // comments).
CFLAGS = -I. -I../ -g  $(PROFILE) $(EXTRA) -g -W -Wall -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wstrict-prototypes -Wmissing-prototypes -Wredundant-decls -Wnested-externs -Wtraditional -Wconversion -Wwrite-strings

LIBS            = -lm 

OFILES          = ma_fitness.o\
		ran1.o\
		poidev.o\
		gammln.o\
		gamdev.o

OBJS		= ma_fitness

${OBJS}: ${OFILES} 
	${CC} ${CFLAGS} -o $@ ${OFILES} ${LIBS}

ma_fitness.o: ma_fitness.c

clean:
	rm -f core ${OBJS} *.o
