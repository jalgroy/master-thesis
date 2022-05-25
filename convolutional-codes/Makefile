CXX := g++
CXXFLAGS := -Wall -g -std=c++17 -pthread -O3 -Wno-sign-compare

srcfiles := channels.cpp
outdir := bin

simulate: $(srcfiles) simulate.cpp
	    $(CXX) $(CXXFLAGS) -o $(outdir)/simulate $(srcfiles) simulate.cpp

profiling: $(srcfiles) simulate.cpp
	    $(CXX) $(CXXFLAGS) -pg -o $(outdir)/simulate_prof $(srcfiles) simulate.cpp

clean:
	    rm -f $(outdir)/$(appname)
