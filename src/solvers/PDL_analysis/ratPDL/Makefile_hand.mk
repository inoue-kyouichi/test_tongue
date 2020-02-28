
##################################################################################
#
# BDIM-cpp

# Copyright (c) 2019-Mechanical and Bioengineering Systems Lab.,
#                    Department of Mechanical Science and Bioengineering,
#                    Graduate School of Engineering Science,
#                    Osaka University.
# All rights reserved.
#
###################################################################################

include ../../../make_setting.mk

TARGET = libPDL_rat.a
TARGET_DIR = ../../lib

CXXSRCS = $(wildcard *.cpp)

SRCS = $(CXXSRCS)

.SUFFIXES: .o .cpp
CXXOBJS = $(CXXSRCS:.cpp=.o)
OBJS  = $(CXXOBJS)
CXXFLAGS +=  -I../. -I../../base -I../../linearSolver -I../../FEM

$(TARGET):$(OBJS)
	-mkdir -p $(TARGET_DIR)
	$(AR) $(TARGET) $(OBJS)
	$(RANLIB) $(TARGET)
	mv $(TARGET) $(TARGET_DIR)

.cpp.o:
	$(CXX) $(CXXFLAGS) $(UDEF_INC_PATH) -o $@ -c $<

clean:
	$(RM) $(OBJS) $(TARGET)

depend: $(CXXOBJS:.o=.cpp)
	@ rm -rf depend.inc
	@ for i in $^; do\
		$(CXX) $(CXXFLAGS) $(UDEF_INC_PATH) -MM $$i >> depend.inc;\
	done

-include depend.inc