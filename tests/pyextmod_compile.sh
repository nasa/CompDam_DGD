#!/usr/bin/env bash

#
# NOTE: must run with -i argument to allow bash to source ~/.bashrc
#

source ~/.bashrc

source activate pyextmod
cd ../pyextmod
make clean
make $MAKE_F90
