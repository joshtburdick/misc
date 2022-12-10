# FIXME: this doesn't track dependencies correctly

# this should work on Linux
# CXXFLAGS = -O3 -Wall -g

# on OSX, with libboost installed via brew, in /opt/homebrew/include/
CXXFLAGS = -O3 -Wall -g -I/opt/homebrew/include/

OBJS = CliqueList.o SubsetIterator.o MaxCliqueCounter.o \
  CliqueCoverageCounter.o utils.o \
  CliquesByEdge.o ZeroableEdgeCliqueCounter.o

all: coverage zeroable_edges

coverage: $(OBJS) coverage.cc
	g++ $(CXXFLAGS) -o coverage $(OBJS) coverage.cc
zeroable_edges: $(OBJS) zeroable_edges.cc
	g++ $(CXXFLAGS) -o zeroable_edges $(OBJS) zeroable_edges.cc


# $(OBJS)

clean:
	rm coverage *.o

