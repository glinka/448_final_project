toCompile = voterModel.cc

CXX = g++

CXXFLAGS = -g -Wall -std=c++0x

all: voterModel

voterModel: $(toCompile)
	$(CXX) $(CXXFLAGS) -o $@ $^

clean:
	$(RM) *.o 

smallSlowTest:
	voterModel 500 4 2 0.5 0.5 0.5 random 500000 10 voterStats

smallFastTest:
	voterModel 500 4 2 0.5 0.5 0.7 random 500000 10 voterStats

bigTest:
	voterModel 10000 4 2 0.5 0.5 0.5 random 1000000 1000 voterStats

debug:
	gdb --args voterModel 100 4 2 0.5 0.5 0.5 random 1000 10 voterStats

