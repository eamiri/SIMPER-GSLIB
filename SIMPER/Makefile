
appname := SIMPER.exe

CXX      := g++ 
CXXFLAGS := -Wno-deprecated -O3 -std=c++11

srcfiles := $(shell find . -name "*.cpp")
objects  := $(patsubst %.cpp, %.o, $(srcfiles))

# include netcdf
#CXXFLAGS += -Dnetcdf                # if netcdf is installed use "CXXFLAGS += -Dnetcdf",                    else "CXXFLAGS += "
#LDLIBS   := -L/usr/local -lnetcdf   # if netcdf is installed give first path "-L<PATH>"and then "-lnetcdf", else "LDLIBS   := "
CXXFLAGS += 
LDLIBS   := 

# if you use a UNIX system, uncomment the LDFLAGS line below
# this is to allow for use a 1Gb, see http://linuxtoosx.blogspot.ca/2010/10/stack-overflow-increasing-stack-limit.html
# LDFLAGS  := -Wl,-stack_size,0x80000000,-stack_addr,0xf0000000

all: $(appname)

$(appname): $(objects)
	$(CXX) $(CXXFLAGS) -o $(appname) $(objects) $(LDLIBS) $(LDFLAGS) 

depend: .depend

.depend: $(srcfiles)
	rm -f ./.depend
	$(CXX) $(CXXFLAGS) -MM $^>>./.depend;

clean:
	rm -f $(objects)

dist-clean: clean
	rm -f *~ .depend
	
all: clean $(binaries)

include .depend

