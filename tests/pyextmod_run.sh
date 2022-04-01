#!/usr/bin/env bash

#
# NOTE: must run with -i argument to allow bash to source ~/.bashrc
#

source ~/.bashrc

source activate compdam
cd ../pyextmod
python verify_debug.py $1 ../tests/testOutput/$2
