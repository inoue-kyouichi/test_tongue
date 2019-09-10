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

.PHONY: depend clean all

all:
	(cd base; make -f Makefile_hand.mk)
	(cd FEM; make -f Makefile_hand.mk)
	(cd RBD; make -f Makefile_hand.mk)
	(cd linearSolver; make -f Makefile_hand.mk)
	(make -f Makefile_PDL.mk)

clean:
	(cd base; make -f Makefile_hand.mk clean)
	(cd FEM; make -f Makefile_hand.mk clean)
	(cd RBD; make -f Makefile_hand.mk clean)
	(cd linearSolver; make -f Makefile_hand.mk clean)
	(make -f Makefile_PDL.mk clean)

depend:
	(cd base; make -f Makefile_hand.mk depend)
	(cd FEM; make -f Makefile_hand.mk depend)
	(cd RBD; make -f Makefile_hand.mk depend)
	(cd linearSolver; make -f Makefile_hand.mk depend)
	(make -f Makefile_PDL.mk depend)