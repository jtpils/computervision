all: greetings cfiles main bye

greetings:
	@echo "::::::::::: Computer vision C++ framework :::::::::::\n"

CCX = g++
CC = gcc

IDIR = ../include
CFLAGS = -I$(IDIR)

ODIR = ../build

LIBS =  -lGL -lGLU -lglfw3 -lX11 -lXxf86vm -lXrandr -lpthread -lXi -ldl \
        -L../lib -lANN -l:libopencv_core.so.4.1 -l:libopencv_imgproc.so.4.1 \
				-l:libopencv_highgui.so.4.1 -l:libopencv_ml.so.4.1 -l:libopencv_video.so.4.1 \
 				-l:libopencv_imgcodecs.so.4.1 -l:libopencv_features2d.so.4.1 \
				-l:libopencv_calib3d.so.4.1 -l:libopencv_objdetect.so.4.1 \
 				-l:libopencv_stitching.so.4.1 \

LDFLAGS = "-Wl,-rpath,../lib/"

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
	$(CCX) -o $@ $^ $(CFLAGS) $(LIBS) $(LDFLAGS)

bye:
	@echo "\nCompilation completed!\n"

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~
