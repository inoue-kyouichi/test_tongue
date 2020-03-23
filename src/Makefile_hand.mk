##################################################################################
#
# BDIM-cpp
#
# Copyright (c) 2019-  Mechanical and Bioengineering Systems Lab.,
#                      Department of Mechanical Science and Bioengineering,
#                      Graduate School of Engineering Science,
#                      Osaka University.
# All rights reserved.
#
###################################################################################

include make_setting.mk
.PHONY: depend clean all utilities install

all:
	(cd solvers; make -f Makefile_hand.mk)
	(cd applications; make -f Makefile_hand.mk)

install:
	cp ../bin/* $(INSTALL_BIN)/

utilities:
	(cd utilities; make -f Makefile_hand.mk utilities)

clean:
	(cd solvers; make -f Makefile_hand.mk clean)
	(cd applications; make -f Makefile_hand.mk clean)
	(cd utilities; make -f Makefile_hand.mk clean)

depend:
	(cd solvers; make -f Makefile_hand.mk depend)
	(cd applications; make -f Makefile_hand.mk depend)
	(cd utilities; make -f Makefile_hand.mk depend)