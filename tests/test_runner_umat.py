#
# Unittest code to run tests on single element models with DGD with Abaqus/Standard and the VUMAT wrapper
#

import os
import shutil
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


    def test_C3D8R_UMAT_matrixTension(self):
        """ Matrix tension, Abaqus/Standard """
        self.runTest("test_C3D8R_UMAT_matrixTension")


    def test_C3D8R_UMAT_simpleShear12(self):
        """ Simple shear in the 1-2 plane, Abaqus/Standard """
        self.runTest("test_C3D8R_UMAT_simpleShear12")



if __name__ == "__main__":
    av.runTests(relPathToUserSub='../for/UMAT', double=True)
