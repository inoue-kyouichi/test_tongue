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
	(cd main_rigidBodyInteraction; make -f Makefile_PDL.mk)

clean:
	(cd main_rigidBodyInteraction; make -f Makefile_PDL.mk clean)

depend:
	(cd main_rigidBodyInteraction; make -f Makefile_PDL.mk depend)