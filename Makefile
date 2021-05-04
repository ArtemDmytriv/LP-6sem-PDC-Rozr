SRCDIR=source/
OBJDIR=obj/
INCDIR=include/
RESDIR=result/

CXX=mpic++
CXXFLAGS+=-std=c++17 -I$(INCDIR)

SOURCES  = $(wildcard $(SRCDIR)*.cpp)
_OBJECTS = $(patsubst $(SRCDIR)%.cpp, %.o, $(SOURCES))
OBJECTS	= $(addprefix $(OBJDIR), $(_OBJECTS))

PROJ=rozr

all: $(PROJ)

$(PROJ): $(OBJECTS) $(RESDIR)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJECTS)

$(OBJECTS): $(OBJDIR)%.o: $(SRCDIR)%.cpp $(OBJDIR)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

$(OBJDIR):
	mkdir $@

$(RESDIR):
	mkdir $@

.PHONY: clean

clean:
	rm -f $(OBJDIR)*.o $(PROJ) result/*.txt *.txt