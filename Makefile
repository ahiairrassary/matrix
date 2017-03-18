EXEC = m

SRC = $(wildcard *.cpp)
OBJ = $(SRC:.cpp=.o)

CXX = g++

CXXFLAGS = -Wall -Wextra -Werror -Wfatal-errors -pedantic -pedantic-errors -Wshadow -Wconversion
LDFLAGS = -lpthread -lblas -llapack


all: $(EXEC)

$(EXEC): $(OBJ)
	$(CXX) -o $(EXEC) $(OBJ) $(LDFLAGS)

%.o: %.cpp
	$(CXX) -o $@ -c $< $(CXXFLAGS) -O3 -std=c++14


.PHONY: clean mrproper

clean:
	rm -rf $(OBJ)

mrproper: clean
	rm -rf $(EXEC)

