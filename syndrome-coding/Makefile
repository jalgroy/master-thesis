CXX := g++
CXXFLAGS := -Wall -std=c++17 -pthread -O3 -Wno-sign-compare

srcfiles := channels.cpp simulate_syndrome.cpp
outdir := bin

appname := simulate
$(appname): $(srcfiles)
	    $(CXX) $(CXXFLAGS) -o $(outdir)/$(appname) $(srcfiles)

clean:
	    rm -f $(outdir)/$(appname)
