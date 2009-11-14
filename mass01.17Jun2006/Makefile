# Project: Mass01
# Makefile created by Dev-C++ 4.9.9.2

CPP  = g++
CC   = gcc
OBJ  = main.o
LINKOBJ  = main.o
LIBS =  -lSDLmain -lSDL  -s 
BIN  = mass01
CXXFLAGS = $(CXXINCS)   -fexpensive-optimizations -O2
CFLAGS = $(INCS) -Dmain=SDL_main   -fexpensive-optimizations -O2
RM = rm -f

.PHONY: all all-before all-after clean clean-custom

all: all-before $(BIN) all-after


clean: clean-custom
	${RM} $(OBJ) $(BIN)

$(BIN): $(OBJ)
	$(CPP) $(LINKOBJ) -o "mass01" $(LIBS)

main.o: main.cc
	$(CPP) -c main.cc -o main.o $(CXXFLAGS)

