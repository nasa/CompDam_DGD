#
# Unittest code to run tests on single element models with DGD with Abaqus/Standard and the VUMAT wrapper
#

import os
import shutil
import re
import abaverify as av

# Helper for dealing with .props files
def copyMatProps():
    # If testOutput doesn't exist, create it
    testOutputPath = os.path.join(os.getcwd(), 'testOutput')
    if not os.path.isdir(testOutputPath):
        os.makedirs(testOutputPath)

    # Put a copy of the properties file in the testOutput directory
    propsFiles = [x for x in os.listdir(os.getcwd()) if x.endswith('.props')]
    for propsFile in propsFiles:
        shutil.copyfile(os.path.join(os.getcwd(), propsFile), os.path.join(os.getcwd(),'testOutput', propsFile))


def modifyParametersFile(jobName='CompDam', **kwargs):
    '''
    For modifying the parameters file
    Input dictionary should have key, value pairs that correspond to entries in CompDam.parameters
    '''

    # Copy/modify parameters file
    with open(os.path.join(os.getcwd(), 'CompDam.parameters'), 'r') as f:
        data = f.read()

    for key, value in kwargs.items():
        data = re.sub(key + r' ?= ?[-0-9\.d(TRUE)(FALSE)]+', key + ' = ' + value, data)

    # Write to testOutput directory
    with open(os.path.join(os.getcwd(), 'testOutput', jobName + '.parameters'), 'w') as f:
        f.write(data)


class SingleElementTests(av.TestCase):
    """
    Single element models to tests the DGD code base with the VUMAT wrapper
    """

    # Class-wide methods
    @classmethod
    def setUpClass(cls):
        copyMatProps()

    # -----------------------------------------------------------------------------------------
    # Test methods
    def test_C3D8R_UMAT_elastic_matrixTension(self):
        """ Elastic matrix tension, Abaqus/Standard """
        self.runTest("test_C3D8R_UMAT_elastic_matrixTension")


    def test_C3D8R_UMAT_elastic_simpleShear12(self):
        """ Elastic simple shear, Abaqus/Standard """
        self.runTest("test_C3D8R_UMAT_elastic_simpleShear12")


    def test_C3D8R_UMAT_matrixTension(self):
        """ Matrix tension, Abaqus/Standard """
        self.runTest("test_C3D8R_UMAT_matrixTension")


    def test_C3D8R_UMAT_simpleShear12(self):
        """ Simple shear in the 1-2 plane, Abaqus/Standard """
        self.runTest("test_C3D8R_UMAT_simpleShear12")


    def test_C3D8R_UMAT_nonlinearShear12(self):
        """ Nonlinear shear model, loading and unloading in 1-2 plane, Abaqus/Standard """
        self.runTest("test_C3D8R_UMAT_nonlinearShear12")


    def test_C3D8R_UMAT_nonlinearShear13(self):
        """ Nonlinear shear model, loading and unloading in 1-3 plane, Abaqus/Standard """
        modifyParametersFile(alpha_search = '.FALSE.')
        self.runTest("test_C3D8R_UMAT_nonlinearShear13")


    def test_C3D8R_UMAT_fiberCompression_FKT_12(self):
        """ Fiber kinking model, Abaqus/Standard """
        self.runTest("test_C3D8R_UMAT_fiberCompression_FKT_12")



if __name__ == "__main__":
    av.runTests(relPathToUserSub='../for/UMAT', double=True)
