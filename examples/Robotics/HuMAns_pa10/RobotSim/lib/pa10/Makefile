CC = gcc
CFLAGS = 
LD = ld
LDFLAGS = -i

SCIDIR = /usr/lib/scilab-3.0


INCLUDES = -I${SCIDIR}/routines
LIB = robotmodel.a

LIB_SRC = Inertia.c NLEffects.c Friction.c Tags.c ScilabInterface.c

all : ${LIB}

${LIB} : ${LIB_SRC}
	${CC} -c ${CFLAGS} ${LIB_SRC} ${INCLUDES}
	${LD} ${LDFLAGS} ${LIB_SRC:.c=.o} -o $@ 

clean :
	rm -f ${LIB} ${LIB_SRC:.c=.o} *~
