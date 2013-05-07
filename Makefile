toCompile = voterModel.cc

CXX = g++

CXXFLAGS = -g -Wall -std=c++0x

all: voterModel

voterModel: $(toCompile)
	$(CXX) $(CXXFLAGS) -o $@ $^

clean:
	$(RM) *.o 

test:
	voterModel 100 4 2 0.5 0.5 0.5 random 1000 10 voterStats

debug:
	gdb --args voterModel 100 4 2 0.5 0.5 0.5 random 1000 10 voterStats

