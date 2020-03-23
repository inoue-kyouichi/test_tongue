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

utilities:
	(cd calcSurfaceCurvature; make -f Makefile.mk)

clean:
	(cd calcSurfaceCurvature; make -f Makefile.mk clean)

depend:
	(cd calcSurfaceCurvature; make -f Makefile.mk depend)
