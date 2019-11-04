all: greetings cfiles main bye

greetings:
	@echo "::::::::::: Computer vision C++ framework :::::::::::\n"

CCX = g++
CC = gcc

IDIR = ../include
CFLAGS = -I$(IDIR)

ODIR = ../build

LIBS =  -lGL -lGLU -lglfw3 -lX11 -lXxf86vm -lXrandr -lpthread -lXi -ldl \
        -L../lib -lANN


_DEPS = reader.h gpa.h viz.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = main.o viz.o glad.o gpa.o reader.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

$(ODIR)/%.o: %.cpp $(DEPS)
	$(CCX) -c -o $@ $< $(CFLAGS)

cfiles:
	@echo "Compiling C and C++ files...\n"
	$(CC) -c glad.c -Wall -o $(ODIR)/glad.o

main: $(OBJ)
	@echo "\nLinking object files...\n"
	$(CCX) -o $@ $^ $(CFLAGS) $(LIBS) 

bye:
	@echo "\nCompilation completed!\n"

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~
