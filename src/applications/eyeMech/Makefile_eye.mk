##################################################################################
#
# FEM SOLID library

# Copyright (c) 2020   Biomechanics Lab.,
#                      Department of Mechanical Science and Bioengineering,
#                      Graduate School of Engineering Science,
#                      Osaka University.
# All rights reserved.
#
###################################################################################

#include ./make_setting.mk
include ../../make_setting.mk

TARGET = eyeMechAnalysis

SRCS = $(wildcard *.cpp)

.SUFFIXES: .o .cpp
CXXOBJS = $(SRCS:.cpp=.o)
OBJS  = $(CXXOBJS)

$(TARGET):$(OBJS)
	$(CXX) $(CXXFLAGS) $(UDEF_INC_PATH) -o $(TARGET) $(OBJS) \
	$(UDEF_LIB_PATH) $(LIBS) $(UDEF_LIBS) $(UDEF_LIB_PATH_SPEC) $(LDFLAGS)
	-mkdir -p $(TARGET_DIR)/bin
	mv $(TARGET) $(TARGET_DIR)/bin

.cpp.o:
	$(CXX) $(CXXFLAGS) $(UDEF_OPT) $(UDEF_INC_PATH) -c $<

clean:
	$(RM) $(OBJS) $(TARGET)

depend: $(OBJS:.o=.cpp)
	@ rm -rf depend.inc
	@ for i in $^; do\
		$(CXX) $(CXXFLAGS) $(UDEF_INC_PATH) -MM $$i >> depend.inc;\
	done

-include depend.inc

