SRCS=vote.cc votingModel.cc votingModelCPI.cc fitCurves.cc calcGraphProps.cc
OBJECTS=$(SRCS:.cc=.o)
GE_SRCS=graph_embedding_main.cc votingModel.cc util_fns.cc eigen_solvers.cc calcGraphProps.cc
GE_OBJECTS=$(GE_SRCS:.cc=.o)

CXX = g++

CXXFLAGS = -g -Wall -std=c++0x -O3 -I/home/alexander/local/eigen/

all: graph_embedding vote

%.o: %.c
	$(CXX) $(CXXFLAGS) -c $<

vote: $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^

graph_embedding: $(GE_OBJECTS)
	$(CXX) -o $@ $^  $(CXXFLAGS)

depend: .depend

.depend: $(SRCS) $(GE_SRCS)
	rm -f ./.depend
	$(CXX) -MM -MT $^ $(CXXFLAGS) > ./.depend

clean:
	$(RM) *.o 

include .depend
