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
	(cd PDL_tooth_interactionAnalysis; make -f Makefile_PDL.mk)
	(cd rat_analysis; make -f Makefile_PDL.mk)
	(cd minimumComplianceProblem; make -f Makefile_Design.mk)

clean:
	(cd PDL_tooth_interactionAnalysis; make -f Makefile_PDL.mk clean)
	(cd rat_analysis; make -f Makefile_PDL.mk clean)
	(cd minimumComplianceProblem; make -f Makefile_Design.mk clean)

depend:
	(cd PDL_tooth_interactionAnalysis; make -f Makefile_PDL.mk depend)
	(cd rat_analysis; make -f Makefile_PDL.mk depend)
	(cd minimumComplianceProblem; make -f Makefile_Design.mk depend)
