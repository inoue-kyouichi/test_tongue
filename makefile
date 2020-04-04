##################################################################################
#
# BDIM Library

# Copyright (c) 2016-9 Mechanical and Bioengineering Systems Lab.,
#                      Department of Mechanical Science and Bioengineering,
#                      Graduate School of Engineering Science,
#                      Osaka University.
# All rights reserved.
#
###################################################################################

.PHONY: depend clean

all:
	cd src && make -f Makefile_hand.mk

test:
	cd src && make -f Makefile_hand.mk test

install:
	cd src && make -f Makefile_hand.mk install

depend:
	cd src && make -f Makefile_hand.mk depend

clean:
	cd src && make -f Makefile_hand.mk clean
