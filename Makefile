toCompile = votingModel.o vote.o

CXX = g++

CXXFLAGS = -g -Wall -std=c++0x #-O2

all: vote

%.o: %.c
	$(CXX) $(CXXFLAGS) -c $<

vote: votingModel.o vote.o
	$(CXX) $(CXXFLAGS) -o $@ $^

clean:
	$(RM) *.o 

debug:
	gdb --args vote random 1000 100000 100 0.5 2 0.5 0.5

