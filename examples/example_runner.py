#
# Code to run example problems
#

import os
import shutil
import abaverify as av
import re


def modifyParametersFile(jobName='CompDam', **kwargs):
    '''
    For modifying the parameters file
    Input dictionary should have key, value pairs that correspond to entries in CompDam.parameters
    '''

    # Copy/modify parameters file
    with open(os.path.join(os.getcwd(), 'CompDam.parameters'), 'r') as f:
        data = f.read()

    for key, value in kwargs.items():
        data = re.sub(key + r' ?= ?[-0-9\.d]*', key + ' = ' + value, data)

    # Write to testOutput directory
    with open(os.path.join(os.getcwd(), 'testOutput', jobName + '.parameters'), 'w') as f:
        f.write(data)


class DCB_fatigue(av.TestCase):
    """
    Demonstrates the cohesive fatigue model with double cantilever beam analyses under displacement control
    """

    # Specify meta class
    __metaclass__ = av.ParametricMetaClass

    # Refers to the template input file name
    baseName = "test_DCB_fatigue"

    # disp is the maximum applied displacement during each fatigue cycle
    parameters = {'disp': [1.48, 1.70, 1.92, 2.25]}

    # staticLoad is the maximum reaction force in the initial static loading step
    expectedpy_parameters = {'staticLoad': [1.813, 2.083, 2.352, 2.754]}

    @classmethod
    def setUpClass(cls):
        modifyParametersFile(
                             fatigue_step = '2',
                             fatigue_R_ratio = '0.1d0',
                             cycles_per_increment_init = '1.d-5',
                             cycles_per_increment_mod = '0.1d0',
                             fatigue_damage_min_threshold = '5.d-5',
                             fatigue_damage_max_threshold = '1.d-3'
                             )


if __name__ == "__main__":
    av.runTests(relPathToUserSub='../for/CompDam_DGD', double=True)
