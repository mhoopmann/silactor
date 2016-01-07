#Set these variables if needed
C = gcc
CC = g++
FLAGS = -O3 -static -D_NOSQLITE -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -DGCC

#Paths to supporting software
MSTOOLKITPATH = ../MSToolkit
KRONIKPATH = ../kronik

#Do not touch these variables
LIBPATH = -L$(MSTOOLKITPATH) -L$(KRONIKPATH)
LIBS = -lmstoolkitlite -lkronik
INCLUDE = -I$(MSTOOLKITPATH)/include -I$(KRONIKPATH) -I./MascotParser


#Do not touch these variables
SILACTOR = CPeptideDatabase.o CSILACtor.o
MASCOTPARSER = MascotParser.o


#Make statements
silactor : $(MASCOTPARSER) $(SILACTOR) SILACtor.cpp
	$(CC) $(FLAGS) $(INCLUDE) $(MASCOTPARSER) $(SILACTOR) SILACtor.cpp $(LIBPATH) $(LIBS) -o silactor

clean:
	rm *.o silactor


#MascotParser objects
MascotParser.o : MascotParser/MascotParser.cpp
	$(CC) $(FLAGS) $(INCLUDE) MascotParser/MascotParser.cpp -c


#SILACtor objects
CPeptideDatabase.o : CPeptideDatabase.cpp
	$(CC) $(FLAGS) $(INCLUDE) CPeptideDatabase.cpp -c

CSILACtor.o : CSILACtor.cpp
	$(CC) $(FLAGS) $(INCLUDE) CSILACtor.cpp -c

