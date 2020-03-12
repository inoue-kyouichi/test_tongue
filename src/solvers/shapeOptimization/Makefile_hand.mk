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
	(cd minimumComplianceProblem; make -f Makefile_hand.mk)

clean:
	(cd minimumComplianceProblem; make -f Makefile_hand.mk clean)

depend:
	(cd minimumComplianceProblem; make -f Makefile_hand.mk depend)