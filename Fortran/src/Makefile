#GCAF MAKEFILE

PROGRAM_NAME = gcaf
COMPILER = gfortran
BUILDER = $(COMPILER)
OPTION= -fopenmp -O3 -g -fcheck=all #-Wall -Wunused-dummy-argument -Wno-tabs
COMPILE = $(COMPILER) $(OPTION) -c

f90sources= main.f90 step.f90 initial.f90 functions.f90 record.f90 read.f90 \
	collisions.f90 field.f90 shared.f90
f90objs= $(f90sources:.f90=.o)
f90mods= shared.mod

$(PROGRAM_NAME): $(f90objs)
	$(COMPILER) $(OPTION) $(f90objs) -o $(PROGRAM_NAME)

shared.mod: shared.o shared.f90
	$(COMPILE) shared.f90

shared.o: shared.f90
	$(COMPILE) shared.f90

collisions.o: shared.f90 collisions.f90
	$(COMPILE) collisions.f90

field.o: shared.f90 field.f90
	$(COMPILE) field.f90

read.o: shared.f90 read.f90
	$(COMPILE) read.f90

record.o: shared.f90 record.f90
	$(COMPILE) record.f90

functions.o: shared.f90 functions.f90
	$(COMPILE) functions.f90

initial.o: shared.f90 initial.f90
	$(COMPILE) initial.f90

step.o: shared.f90 step.f90
	$(COMPILE) step.f90

main.o: shared.mod main.f90
	$(COMPILE) main.f90


.PHONY: run clean install export

clean:
	rm -f $(f90objs) $(f90mods) $(PROGRAM_NAME)

install:
	cp gcaf ~/.local/bin

export:
	cp Makefile ~/github/gcaf/src/
	cp -r *.f90 ~/github/gcaf/src/
