# Makefile to compile all cpp and header files in current directory

# Compiler to use
CXX = g++

# Flags to pass to the compiler
CXXFLAGS = -std=c++14 -Wall -I./lib -O3 -march=core-avx2 -pthread -std=c++14


# Find all cpp and header files in current directory
SRCS := $(wildcard *.cpp)
HDRS := $(wildcard *.h)

# Object files to generate
OBJS := $(SRCS:.cpp=.o)

# Name of the final executable
EXEC = cover_tree

all: $(EXEC)

$(EXEC): $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) -o $(EXEC)

%.o: %.cpp $(HDRS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(EXEC)

