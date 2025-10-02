CC = g++
FLAGS = -std=c++17

OUT = Ex05

SOURCES = Ex05.cpp helper.cpp eigen.cpp
OBJECTS = $(SOURCES:.cpp=.o)

all: $(OUT)

$(OUT): $(OBJECTS)
	$(CC) $(OBJECTS) -o $(OUT) $(FLAGS)

%.o: %.cpp
	$(CC) -c $< -o $@ $(FLAGS)

.PHONY: clean
clean:
	rm -f $(OBJECTS) $(OUT)