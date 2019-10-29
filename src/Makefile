IDIR = ../include
CC=g++
CFLAGS=-I$(IDIR)

ODIR= ../build

LIBS= -L../lib -lANN

_DEPS = reader.h gpa.h 
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = reader.o gpa.o main.o ann/ANN.o ann/kd_tree.o ann/perf.o ann/kd_util.o ann/kd_split.o ann/kd_search.o ann/brute.o ann/kd_dump.o ann/bd_tree.o ann/bd_search.o ann/bd_fix_rad_search.o ann/kd_pr_search.o ann/kd_fix_rad_search.o ann/bd_pr_search.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

$(ODIR)/%.o: %.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

main: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~ 
