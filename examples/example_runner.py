#
# Code to run example problems
#

import os
import shutil
import abaverify as av


class UNC0_FKT(av.TestCase):
    """
    Demonstrates the fiber kinking theory model for unnotched compression
    """

    # Specify meta class
    __metaclass__ = av.ParametricMetaClass

    # Refers to the template input file name
    baseName = "test_UNC0_C3D8R_FKT"

    # Use python script instead of input file
    pythonScriptForModel = True

    # Range of parameters to test; all combinations are tested
    parameters = {'meshSize': [0.05, 0.1, 0.2]}




if __name__ == "__main__":
    av.runTests(relPathToUserSub='../for/CompDam_DGD', double=True)
