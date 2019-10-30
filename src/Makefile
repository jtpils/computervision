IDIR = ../include
CC=g++
CFLAGS=-I$(IDIR)

ODIR= ../build

LIBS= -L../lib -lANN

_DEPS = reader.h gpa.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = reader.o gpa.o main.o 
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

$(ODIR)/%.o: %.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

main: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~
