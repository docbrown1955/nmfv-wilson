# Part depending on the generated library

NAME    = nmfw_wilson
NAMELIB = libnmfv

INCPATH = -Iinclude 
LIBPATH = 
LIBS    = -looptools-quad -lgfortran -lquadmath  -lgsl -lgslcblas

# Part common to all libraries generated by MARTY

CXX     = g++
MATH_OPTI = \
		-fno-math-errno\
		-ffinite-math-only\
		-fno-rounding-math\
		-fno-signaling-nans\
		-fcx-limited-range\
		-fexcess-precision=fast\
		-fno-signed-zeros\
		-fno-trapping-math\
		-fassociative-math\
		-freciprocal-math

DEFAULTFLAGS = -std=c++17 -Wall -Wextra -Wpedantic -Wno-deprecated-declarations -fPIC
OPTIFLAGS = # $(MATH_OPTI)
QUADFLAGS = # -DQUAD=1 -DQUADSIZE=16
CXXFLAGS  = $(DEFAULTFLAGS) $(OPTIFLAGS) $(QUADFLAGS)
LINKFLAGS = $(DEFAULTFLAGS) $(QUADFLAGS)

SRCDIR = src
INCDIR = include
OBJDIR = obj
SCRIPTDIR = script
SOBJDIR = script/obj
BINDIR = bin
LIBDIR = lib

SRC 	 = $(wildcard $(SRCDIR)/*.cpp)
HEADERS  = $(wildcard $(INCDIR)/*.h)
OBJ_INIT = $(SRC:$(SRCDIR)/%.cpp=$(OBJDIR)/%.o)
OBJ 	 = $(filter-out $(OBJDIR)/$(NAME)_pylink.o, $(OBJ_INIT))
SCRIPTS  = $(wildcard $(SCRIPTDIR)/*.cpp)
BINARIES = $(SCRIPTS:$(SCRIPTDIR)/%.cpp=%.x)
SOBJ = $(SCRIPTS:$(SCRIPTDIR)/%.cpp=$(SOBJDIR)/%.o)

all: scripts lib

scripts: $(BINARIES)
lib: $(NAMELIB).so $(NAMELIB).a

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@ $(INCPATH)

$(SOBJDIR)/%.o: $(SCRIPTDIR)/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@ $(INCPATH)

%.x: $(SOBJDIR)/%.o $(OBJ)
	$(CXX) $(LINKFLAGS) -o $(BINDIR)/$@ $< $(OBJ) $(INCPATH) $(LIBPATH) $(LIBS)

$(NAMELIB).so: $(OBJ)
	$(CXX) $(CXXFLAGS) -shared -o $(LIBDIR)/$@ $(OBJ) $(LIBPATH) $(LIBS)
$(NAMELIB).a: $(OBJ)
	ar rcs $(LIBDIR)/$@ $(OBJ)

clean:
	rm $(OBJDIR)/*.o
Clean:
	$(MAKE) clean
	rm $(BINDIR)/*.x