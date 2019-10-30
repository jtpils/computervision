all: greetings cfiles viz bye

greetings:
	@echo ":::::: Computer vision C++ framework ::::::\n"

CCX=g++
CC=gcc

IDIR = ../include
CFLAGS=-I$(IDIR)

ODIR= ../build

LIBS=  -lGL -lGLU -lglfw3 -lX11 -lXxf86vm -lXrandr -lpthread -lXi -ldl -L../lib -lANN

_DEPS = reader.h gpa.h viz.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = viz.o glad.o gpa.o reader.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

$(ODIR)/%.o: %.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

cfiles:
	@echo "Compiling C and C++ files...\n"
	$(CCX) -c glad.c -Wall -o $(ODIR)/glad.o
	@echo "\n"

viz: $(OBJ)
	@echo "Linking object files...\n"
	$(CCX) -o $@ $^ $(CFLAGS) $(LIBS)

bye:
	@echo "Compilation completed!\n"

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~
