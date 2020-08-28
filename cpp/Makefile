# GCAFPP MAKEFILE

PROGRAM_NAME = gcafpp
COMPILER = gcc
BOOSTPATH = -I ~/boost_1_70_0
BUILDER = $(COMPILER)
OPT = -O3 -lm -lstdc++ -fopenmp -std=c++14
BOOSTOPT = -L/stage/lib -lboost_program_options
BOOSTOPT_LOG = -lboost_log -lboost_thread -lpthread -lboost_system -lboost_log_setup
COMPILE = $(COMPILER) $(OPT) -c

cppsources= main.cpp read.cpp record.cpp initial.cpp field.cpp step.cpp functions.cpp collisions.cpp
cppobjs= $(cppsources:.cpp=.o)

$(PROGRAM_NAME): $(cppobjs)
	$(COMPILER) $(cppobjs) -o $(PROGRAM_NAME) $(OPT) $(BOOSTOPT)

functions.o: shared.h functions.cpp
	$(COMPILE) functions.cpp

read.o: shared.h read.cpp
	$(COMPILE) $(BOOSTPATH) read.cpp

record.o: shared.h record.cpp
	$(COMPILE) record.cpp

initial.o: shared.h initial.cpp
	$(COMPILE) initial.cpp

field.o: shared.h field.cpp
	$(COMPILE) field.cpp

collisions.o: shared.h collisions.cpp
	$(COMPILE) collisions.cpp

step.o: shared.h step.cpp
	$(COMPILE) step.cpp

main.o: shared.h define.h main.cpp
	$(COMPILE) main.cpp


.PHONY: clean install export xook

clean:
	rm -f $(PROGRAM_NAME) $(cppobjs)
	rm -f *.plt

install:
	cp $(PROGRAM_NAME) ~/.local/bin

export:
	cp Makefile ~/github/gcafpp/src/
	cp -r *.cpp ~/github/gcafpp/src/
	cp -r *.h ~/github/gcafpp/src/

xook:
	scp *.cpp saavedra@xook.nucleares.unam.mx:/storage/icn/saavedra/src
	scp *.h saavedra@xook.nucleares.unam.mx:/storage/icn/saavedra/src
