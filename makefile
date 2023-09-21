GCC = g++
CFLAGS = -std=c++98 -Wall -I/opt/homebrew/include
LDFLAGS = -L/opt/homebrew/lib
LDLIBS = -lhdf5

all: gas

gas: gas.cpp
	$(GCC) $(CFLAGS) $(LDFLAGS) -o gas gas.cpp $(LDLIBS)

run: gas
	./gas
	python plotting.py
	open animation.mp4

clean:
	rm -rf gas
