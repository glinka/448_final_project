toCompile = vote.o voterModelAlphaInitDist.o

CXX = g++

CXXFLAGS = -g -Wall -O2 -std=c++0x

all: voterModelAlphaInitDist

%.o: %.c
	$(CXX) $(CXXFLAGS) -c $<

voterModelAlphaInitDist: voterModelAlphaInitDist.o vote.o
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

