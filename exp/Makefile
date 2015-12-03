CC = g++

CFLAGS = -g -Wall -std=c++11

LDFLAGS =

SOURCES = generate.cpp

OBJECTS = $(SOURCES:.cpp=.o)

EXECUTABLE = generate

all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

clean:
	rm -rf $(EXECUTABLE) $(OBJECTS)
	rm -rf coloring
	
.cpp.o:
	$(CC) -I$(INCDIR) $(CFLAGS) $< -c -o $@