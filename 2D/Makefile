CC 		= g++
CFLAGS 	= -c
STDFLAGS = -std=c++11
EXNAME 	= BreastCancer_2D

OBJS 	= main.o BreastCancer_2D.o Initialization.o Packing.o MD.o Minimizer.o Force.o VerletList.o LocalPhi.o FileIO.o misc.o 

all: ${OBJS}
	${CC} ${STDFLAGS} -o ${EXNAME} ${OBJS}

main.o: main.cpp
	${CC} ${STDFLAGS} ${CFLAGS} main.cpp

BreastCancer_2D.o: BreastCancer_2D.cpp BreastCancer_2D.h
	${CC} ${STDFLAGS} ${CFLAGS} BreastCancer_2D.cpp

Initialization.o: Initialization.cpp BreastCancer_2D.h
	${CC} ${STDFLAGS} ${CFLAGS} Initialization.cpp

Packing.o: Packing.cpp BreastCancer_2D.h
	${CC} ${STDFLAGS} ${CFLAGS} Packing.cpp

MD.o: MD.cpp BreastCancer_2D.h
	${CC} ${STDFLAGS} ${CFLAGS} MD.cpp

Minimizer.o: Minimizer.cpp BreastCancer_2D.h
	${CC} ${STDFLAGS} ${CFLAGS} Minimizer.cpp

Force.o: Force.cpp BreastCancer_2D.h
	${CC} ${STDFLAGS} ${CFLAGS} Force.cpp

VerletList.o: VerletList.cpp BreastCancer_2D.h
	${CC} ${STDFLAGS} ${CFLAGS} VerletList.cpp

LocalPhi.o: LocalPhi.cpp BreastCancer_2D.h
	${CC} ${STDFLAGS} ${CFLAGS} LocalPhi.cpp

FileIO.o: FileIO.cpp BreastCancer_2D.h
	${CC} ${STDFLAGS} ${CFLAGS} FileIO.cpp

misc.o: misc.cpp BreastCancer_2D.h
	${CC} ${STDFLAGS} ${CFLAGS} misc.cpp
	
clean:
	rm -f ${OBJS} !${EXNAME}
	@echo "all objects cleaned up!"
