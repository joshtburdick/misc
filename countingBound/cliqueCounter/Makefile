CXXFLAGS = -O6 -Wall -g
OBJS = CliqueList.o SubsetIterator.o MaxCliqueCounter.o \
  CliqueCoverageCounter.o

all: coverage

coverage: $(OBJS) coverage.cc
	g++ $(CXXFLAGS) -o coverage $(OBJS) coverage.cc


# $(OBJS)

clean:
	rm coverage *.o

