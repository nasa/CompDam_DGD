from abaqus import *
import abaqusConstants
from abaqusConstants import *
from caeModules import *
import mesh
import interaction
import step

import os
import inspect
import re

# Parameters
modelName = 'Model-1'
partName = 'Part-1'
jobName = os.path.basename(inspect.getfile(inspect.currentframe())).split('.')[0]
length = 14
width = 5
numPlies = 6
plyThickness = 0.183
undamageableLength = 4
wkb = 0.05
meshSize = 0.2
stepDuration = 0.1
featureFlags = 10300

# length = 2*fineMeshSize

density = 1.57e-09 * 10000

tol = 1e-3*min(length, width, plyThickness, meshSize)

# Throw-away debugger
def debug(string):
    if True:   # Use to turn off debugging
        sys.__stderr__.write("DEBUG - " + __name__ + ":  " + string + '\n')


def sortEdgesForBiasSeed(part, edges, center):
    '''
    Arguments:
    part: reference to the part object
    edges: a tuple of edges
    center: a tuple of x, y, z coordinates indicating the center of the fine mesh region

    This function returns the tuples e1 and e2 which corresponds to the arguments needed in seedEdgeByBias()
    where end1Edges = e1 and end2Edges = e2. In the seedEdgeByBias(), the smallest elements are
    positioned near the first vertice on the edge for edges passed in to end1Edges and visa versa.
    '''

    e1 = list()
    e2 = list()

    for e in edges:
        (v1, v2) = e.getVertices()
        (v1Coords, v2Coords) = (p.vertices[v1].pointOn[0], p.vertices[v2].pointOn[0])
        v1Distance = sqrt((v1Coords[0]-center[0])**2 + (v1Coords[1]-center[1])**2 + (v1Coords[2]-center[2])**2)
        v2Distance = sqrt((v2Coords[0]-center[0])**2 + (v2Coords[1]-center[1])**2 + (v2Coords[2]-center[2])**2)

        if v1Distance < v2Distance:
            e1.append(e)
        else:
            e2.append(e)

    return (tuple(e1), tuple(e2))

# Create geometry, sets, boundary conditions
mdb.Model(name=modelName, modelType=STANDARD_EXPLICIT)
s = mdb.models[modelName].ConstrainedSketch(name='__profile__', sheetSize=5*max(length, width))
g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
s.rectangle(point1=(0.0, 0.0), point2=(length, width))

# Model definition
mdb.models[modelName].Part(name=partName, dimensionality=THREE_D, type=DEFORMABLE_BODY)
m = mdb.models[modelName]
p = m.parts[partName]
p.BaseShell(sketch=s)

# Partition the model
p.DatumPlaneByPrincipalPlane(offset=undamageableLength/2, principalPlane=YZPLANE)
p.DatumPlaneByPrincipalPlane(offset=length-undamageableLength/2, principalPlane=YZPLANE)
for plane in p.datums.keys():
    p.PartitionFaceByDatumPlane(faces=p.faces, datumPlane=p.datums[plane])

# Materials
m.Material(name='IM7-8552')
m.materials['IM7-8552'].UserMaterial(mechanicalConstants=( \
    featureFlags, 0, plyThickness,  0,  0,  0,  0,  0, \
    152689.0, 8703.0,   5164.0,   0.32,     0.45,     62.3,     92.30,    0.277,    \
    0.788,    1.634,    199.8,    0.925,      0,          0,         0,         0,   \
    -5.5e-6,  2.58e-5,  4.06e-9,    5.4,     2326.2,   0.2,      133.3,    0.5, \
    1731.,      0,        0,         0,          0,    wkb,     0,   0.3))
m.materials['IM7-8552'].Density(table=((density, ), ))

m.Material(name='IM7-8552-Elastic')
m.materials['IM7-8552-Elastic'].UserMaterial(mechanicalConstants=( \
    0, 0, plyThickness,  0,  0,  0,  0,  0, \
    152689.0, 8703.0,   5164.0,   0.32,     0.45,     62.3,     92.30,    0.277,    \
    0.788,    1.634,    199.8,    0.925,      0,          0,         0,         0,   \
    -5.5e-6,  2.58e-5,  4.06e-9,    5.4,     2326.2,   0.2,      133.3,    0.5, \
    1731.,      0,        0,         0,          0,    wkb,     0,   0.3))
m.materials['IM7-8552-Elastic'].Density(table=((density, ), ))

# Create seeds
p.seedEdgeBySize(edges=p.edges, size=meshSize)

# Mesh controls
p.setMeshControls(technique=STRUCTURED, elemShape=QUAD, regions=p.faces)

# Generate mesh
p.generateMesh()

# Create the mesh part
p.PartFromMesh(copySets=True, name='orphanMeshPart')
p = m.parts['orphanMeshPart']

# Extrude shell layer for facesheet
p.generateMeshByOffset(meshType=SOLID, numLayers=numPlies, region=regionToolset.Region(side1Elements=p.elements),
    shareNodes=True, totalThickness=plyThickness*numPlies, deleteBaseElements=True)

# Sets for material properties
e1 = p.elements.getByBoundingBox(xMax=undamageableLength/2 + tol)
e2 = p.elements.getByBoundingBox(xMin=length-undamageableLength/2 - tol)
p.Set('elastic', elements=e1+e2)
e = p.elements.getByBoundingBox(xMin=undamageableLength/2 - tol, xMax=length-undamageableLength/2 + tol)
p.Set('damageable', elements=e)


p.Set('all', elements=p.elements)
p.Set('allNodes', elements=p.nodes)

# Properties
m.HomogeneousSolidSection(name='Section-damageable', material='IM7-8552', thickness=None)
p.SectionAssignment(region=p.sets['damageable'], sectionName='Section-damageable', thicknessAssignment=FROM_SECTION)

m.HomogeneousSolidSection(name='Section-elastic', material='IM7-8552-Elastic', thickness=None)
p.SectionAssignment(region=p.sets['elastic'], sectionName='Section-elastic', thicknessAssignment=FROM_SECTION)

# Material orientation
csys = p.DatumCsysByThreePoints(coordSysType=CARTESIAN, line1=(1.0, 0.0, 0.0), line2=(0.0, 1.0, 0.0), name='PartCsys', origin=(0.0, 0.0, 0.0))
p.MaterialOrientation(axis=AXIS_3, localCsys=p.datums[csys.id], orientationType=SYSTEM, region=p.sets['all'], stackDirection=STACK_3)

# Sets for boundary conditions
# Pinned
n = p.nodes.getByBoundingBox(xMax=tol, yMax=tol, zMax=tol)
p.Set('BC-pin', nodes=n)
# Clamp left
n = p.nodes.getByBoundingBox(xMax=tol)
leftSide = p.Set('BC-clamp-left', nodes=n)
# Clamp right
n = p.nodes.getByBoundingBox(xMin=length - tol)
rightSide = p.Set('BC-clamp-right', nodes=n)
# Load app
n = p.nodes.getByBoundingBox(xMin=length - tol, yMax=tol, zMax=tol)
loadAppRight = p.Set('BC-loadApp-right', nodes=n)
n = p.nodes.getByBoundingBox(xMax=tol, yMax=tol, zMax=tol)
loadAppLeft = p.Set('BC-loadApp-left', nodes=n)
# Load app followers
p.SetByBoolean(name='BC-loadAppFollowers-right', sets=[rightSide, loadAppRight], operation=DIFFERENCE)
p.SetByBoolean(name='BC-loadAppFollowers-left', sets=[leftSide, loadAppLeft], operation=DIFFERENCE)

# Assign element type
p.Set('all', elements=p.elements)
p.setElementType(elemTypes=(mesh.ElemType(elemCode=C3D8R, elemLibrary=EXPLICIT, secondOrderAccuracy=OFF, hourglassControl=ENHANCED), ), regions=(p.sets['elastic']))
p.setElementType(elemTypes=(mesh.ElemType(elemCode=C3D8R, elemLibrary=EXPLICIT, secondOrderAccuracy=OFF, hourglassControl=ENHANCED), ), regions=(p.sets['damageable']))

# Assembly
instanceName = 'orphanMeshPart-1'
m.rootAssembly.DatumCsysByDefault(CARTESIAN)
m.rootAssembly.Instance(dependent=ON, name=instanceName, part=p)

# Step
m.ExplicitDynamicsStep(name='Step-1', previous='Initial', timePeriod=stepDuration, massScaling=((SEMI_AUTOMATIC, MODEL, THROUGHOUT_STEP, 0.0, 5.0e-07, BELOW_MIN, 1, 0, 0.0, 0.0, 0, None), ), nlgeom=ON)

# Amplitude
m.SmoothStepAmplitude(name='Amp-1', timeSpan=STEP, data=((0.0, 0.0), (stepDuration, 1.0)))

# Boundary conditions
m.DisplacementBC(name='BC-clamp-Left', createStepName='Step-1', region=m.rootAssembly.instances[instanceName].sets['BC-clamp-left'], u1=UNSET, u2=0.0, u3=0.0, amplitude='Amp-1',)
m.DisplacementBC(name='BC-loadApp-Left', createStepName='Step-1', region=m.rootAssembly.instances[instanceName].sets['BC-loadApp-left'], u1=0.08, u2=UNSET, u3=UNSET, amplitude='Amp-1',)
m.Equation(name='LoadApp-left', terms=((-1.0, instanceName + '.BC-loadAppFollowers-left', 1), (1.0, instanceName + '.BC-loadApp-left', 1)))
m.DisplacementBC(name='BC-clamp-Right', createStepName='Step-1', region=m.rootAssembly.instances[instanceName].sets['BC-clamp-right'], u1=UNSET, u2=0.0, u3=0.0, amplitude='Amp-1',)
m.DisplacementBC(name='BC-loadApp-Right', createStepName='Step-1', region=m.rootAssembly.instances[instanceName].sets['BC-loadApp-right'], u1=-0.08, u2=UNSET, u3=UNSET, amplitude='Amp-1',)
m.Equation(name='LoadApp-right', terms=((-1.0, instanceName + '.BC-loadAppFollowers-right', 1), (1.0, instanceName + '.BC-loadApp-right', 1)))

# Field output
m.FieldOutputRequest(createStepName='Step-1', numIntervals=50, name='all', region=m.rootAssembly.instances[instanceName].sets['all'], sectionPoints=DEFAULT, variables=('S', 'LE', 'SDV'))
m.FieldOutputRequest(createStepName='Step-1', numIntervals=50, name='all-nodes', region=m.rootAssembly.instances[instanceName].sets['allNodes'], sectionPoints=DEFAULT, variables=('U', ))

del m.fieldOutputRequests['F-Output-1']


# History output - load app
m.HistoryOutputRequest(createStepName='Step-1', frequency=1, name='LoadApp-left', region=m.rootAssembly.instances[instanceName].sets['BC-loadApp-left'], sectionPoints=DEFAULT, variables=('U1', 'RF1'))
m.HistoryOutputRequest(createStepName='Step-1', frequency=1, name='LoadApp-right', region=m.rootAssembly.instances[instanceName].sets['BC-loadApp-right'], sectionPoints=DEFAULT, variables=('U1', 'RF1'))

# Step 2 for detailed output at failure
m.ExplicitDynamicsStep(name='Step-2', previous='Step-1', timePeriod=stepDuration )
m.boundaryConditions['BC-loadApp-Left'].setValuesInStep(stepName='Step-2', u1=0.01)
m.boundaryConditions['BC-loadApp-Right'].setValuesInStep(stepName='Step-2', u1=-0.01)
m.fieldOutputRequests['all'].setValuesInStep(stepName='Step-2', numIntervals=500)
m.fieldOutputRequests['all-nodes'].setValuesInStep(stepName='Step-2', numIntervals=500)

# Unsupported keywords
m.keywordBlock.synchVersions(storeNodesAndElements=False)
position = m.keywordBlock.sieBlocks.index('*Material, name=IM7-8552-Elastic') + 2
m.keywordBlock.insert(position,
    '*Depvar\n' +
    ' 24,\n' +
    '  1, CDM_d2\n' +
    '  2, CDM_Fb1\n' +
    '  3, CDM_Fb2\n' +
    '  4, CDM_Fb3\n' +
    '  5, CDM_B\n' +
    '  6, CDM_Lc\n' +
    '  7, CDM_FIfT\n' +
    '  8, CDM_d1\n' +
    '  9, CDM_FIm\n' +
    ' 10, CDM_alpha\n' +
    ' 11, CDM_STATUS\n' +
    ' 12, CDM_Plas12\n' +
    ' 13, CDM_Inel12\n' +
    ' 14, CDM_eps12\n' +
    ' 15, CDM_slide1\n' +
    ' 16, CDM_slide2\n' +
    ' 17, CDM_FIfC\n' +
    ' 18, CDM_d1T\n' +
    ' 19, CDM_d1C\n' +
    ' 20, CDM_phi\n' +
    ' 21, CDM_gamma\n' +
    ' 22, CDM_Fm1\n' +
    ' 23, CDM_Fm2\n' +
    ' 24, CDM_Fm3\n' +
    '*Characteristic Length, definition=USER, components=6\n' +
    '*Initial Conditions, Type=Solution\n'+
    ' orphanMeshPart-1.elastic,  0.d0,  0.d0,  0.d0,  0.d0,  0.d0,  0.d0,  0.d0,\n'+
    ' 0.d0,  0.d0,  -999,     1,  0.d0,  0.d0,  0.d0,  0.d0,\n'+
    ' 0.d0,  0.d0,  0.d0,  0.d0, 0.d0,  0.d0,  0.d0,  0.d0,\n'+
    ' 0.d0,')

position = m.keywordBlock.sieBlocks.index('*Material, name=IM7-8552') + 2
m.keywordBlock.insert(position,
    '*Depvar\n' +
    ' 24,\n' +
    '  1, CDM_d2\n' +
    '  2, CDM_Fb1\n' +
    '  3, CDM_Fb2\n' +
    '  4, CDM_Fb3\n' +
    '  5, CDM_B\n' +
    '  6, CDM_Lc\n' +
    '  7, CDM_FIfT\n' +
    '  8, CDM_d1\n' +
    '  9, CDM_FIm\n' +
    ' 10, CDM_alpha\n' +
    ' 11, CDM_STATUS\n' +
    ' 12, CDM_Plas12\n' +
    ' 13, CDM_Inel12\n' +
    ' 14, CDM_eps12\n' +
    ' 15, CDM_slide1\n' +
    ' 16, CDM_slide2\n' +
    ' 17, CDM_FIfC\n' +
    ' 18, CDM_d1T\n' +
    ' 19, CDM_d1C\n' +
    ' 20, CDM_phi\n' +
    ' 21, CDM_gamma\n' +
    ' 22, CDM_Fm1\n' +
    ' 23, CDM_Fm2\n' +
    ' 24, CDM_Fm3\n' +
    '*Characteristic Length, definition=USER, components=6\n' +
    '*Initial Conditions, Type=Solution\n'+
    ' orphanMeshPart-1.damageable,  0.d0,  0.d0,  0.d0,  0.d0,  0.d0,  0.d0,  0.d0,\n'+
    ' 0.d0,  0.d0,  -999,     1,  0.d0,  0.d0,  0.d0,  0.d0,\n'+
    ' 0.d0,  0.d0,  0.d0,  0.d0,  1.0d0,  0.d0,  0.d0,  0.d0,\n'+
    ' 0.d0,')


# Create job
mdb.Job(model=modelName, name=jobName)
mdb.jobs[jobName].writeInput()
