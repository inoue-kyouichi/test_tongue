COMPILER = icpc

#local environment
DEST   = ~/bin
LIBDIR = /home/totani/lib
TARDIR = ./bin
OBJDIR = ./obj
SRCDIR = ./src

CPPFLAGS = -Wall -fast -std=c++0x -g -MMD -MP -fopenmp
LDFLAGS  = -fopenmp
LIBS     =
INCLUDE  = -I./include

#MKL
LDFLAGS  += -mkl

#Eigen
INCLUDE  += -I$(LIBDIR)/eigen-3.3.4

#TEXTPARSER
TXTP      = $(LIBDIR)/TextParser-1.8.5
LDFLAGS  += -L$(TXTP)/lib -lTP
INCLUDE  += -I$(TXTP)/include

#glog
GLOG      = $(LIBDIR)/glog
LDFLAGS += -L$(GLOG)/lib -lglog
INCLUDE += -I$(GLOG)/include

TARGET  = $(TARDIR)/femAnalysis

ifeq "$(strip $(SRCDIR))" ""
  SRCDIR = .
endif
ifeq "$(strip $(OBJDIR))" ""
  OBJDIR = .
endif

SOURCES = $(wildcard $(SRCDIR)/*.cpp)
OBJECTS = $(addprefix $(OBJDIR)/, $(notdir $(SOURCES:.cpp=.o)))
DEPENDS = $(OBJECTS:.o=d)


$(TARGET): $(OBJECTS) $(LIBS)
	@[ -d $(TARDIR) ] || mkdir -p $(TARDIR)
	$(COMPILER) -o $@ $^ $(LDFLAGS)

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	@[ -d $(OBJDIR) ] || mkdir -p $(OBJDIR)
	$(COMPILER) $(CPPFLAGS) $(INCLUDE) -o $@ -c $<

all: clean $(TARGET)

clean:
	rm -f $(OBJECTS) $(DEPENDS) $(TARGET)

install:$(TARGET)
	install -s $(TARGET) $(DEST)

-include $(DEPENDS)
