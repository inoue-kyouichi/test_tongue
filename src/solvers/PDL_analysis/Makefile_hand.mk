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
	(cd PDL; make -f Makefile_hand.mk)
	(cd PDL_TOOTH_interaction; make -f Makefile_hand.mk)

clean:
	(cd PDL; make -f Makefile_hand.mk clean)
	(cd PDL_TOOTH_interaction; make -f Makefile_hand.mk clean)

depend:
	(cd PDL; make -f Makefile_hand.mk depend)
	(cd PDL_TOOTH_interaction; make -f Makefile_hand.mk depend)