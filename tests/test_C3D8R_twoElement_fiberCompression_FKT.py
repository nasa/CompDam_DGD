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
length = 0.2
width = 0.2
thickness = 0.15
elasticElToTotal = 0.5
# fineMeshSize = 0.1
# coarseMeshSize = 0.5
# fineMeshExtent = 0.4
stepDuration = 0.1
featureFlags = 10300

# length = 2*fineMeshSize

density = 1.57e-09 * 1000

tol = 1e-3*min(length, width, thickness)

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
p.DatumPlaneByPrincipalPlane(offset=length-elasticElToTotal*length, principalPlane=YZPLANE)
for plane in p.datums.keys():
    p.PartitionFaceByDatumPlane(faces=p.faces, datumPlane=p.datums[plane])

# Materials
m.Material(name='IM7-8552')
m.materials['IM7-8552'].UserMaterial(mechanicalConstants=( \
    featureFlags, 0, thickness,  0,  0,  0,  0,  0, \
    171420.0, 9080.0,   5290.0,   0.32,     0.52,     62.3,     92.30,    0.277,    \
    0.788,    1.634,    199.8,    0.925,      0,          0,         0,         0,   \
    -5.5e-6,  2.58e-5,  4.412e-10, 5.934,     2326.2,   0.2,      133.3,    0.5, \
    1200.1,      0,        0,         0,          0,  length,         0,    0.3))
m.materials['IM7-8552'].Density(table=((density, ), ))
m.Material(name='IM7-8552-Elastic')
m.materials['IM7-8552-Elastic'].Elastic(table=((171420.0, 9080.0, 9080.0, 0.32, 0.32, 0.52, 5290.0, 5290.0, 2986.8), ), type=ENGINEERING_CONSTANTS)
m.materials['IM7-8552-Elastic'].Density(table=((density, ), ))

# Create seeds
p.seedEdgeByNumber(edges=p.edges, number=1)

# Mesh controls
p.setMeshControls(technique=STRUCTURED, elemShape=QUAD, regions=p.faces)

# Generate mesh
p.generateMesh()

# Create the mesh part
p.PartFromMesh(copySets=True, name='orphanMeshPart')
p = m.parts['orphanMeshPart']

# Extrude shell layer for facesheet
p.generateMeshByOffset(meshType=SOLID, numLayers=1, region=regionToolset.Region(side1Elements=p.elements),
    shareNodes=True, totalThickness=thickness, deleteBaseElements=True)

# Sets for material properties
e = p.elements.getByBoundingBox(xMax=length-elasticElToTotal*length + tol)
p.Set('elastic', elements=e)
e = p.elements.getByBoundingBox(xMin=length-elasticElToTotal*length - tol)
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
p.MaterialOrientation(axis=AXIS_3, localCsys=p.datums[csys.id], orientationType=SYSTEM, region=p.sets['damageable'], stackDirection=STACK_3)
p.MaterialOrientation(axis=AXIS_3, localCsys=p.datums[csys.id], orientationType=SYSTEM, region=p.sets['elastic'], stackDirection=STACK_3, additionalRotationType=ROTATION_ANGLE, angle=2.167)

# Sets for boundary conditions
# Pinned
n = p.nodes.getByBoundingBox(xMax=tol, yMax=tol, zMax=tol)
p.Set('BC-pin', nodes=n)
# Clamp left
n = p.nodes.getByBoundingBox(xMax=tol)
p.Set('BC-clamp-left', nodes=n)
# Clamp right
n = p.nodes.getByBoundingBox(xMin=length - tol)
temp = p.Set('BC-clamp-right', nodes=n)
# Load app
n = p.nodes.getByBoundingBox(xMin=length - tol, yMax=tol, zMax=tol)
loadApp = p.Set('BC-loadApp', nodes=n)
# Load app followers
p.SetByBoolean(name='BC-loadAppFollowers', sets=[temp, loadApp], operation=DIFFERENCE)
# Prevent excessive mesh distortion
midNodes = p.nodes.getByBoundingBox(xMin=length-elasticElToTotal*length - tol, xMax=length-elasticElToTotal*length + tol)
n = midNodes.getByBoundingBox(yMin=width-tol, zMin=thickness-tol)
p.Set('BC-left-upper', nodes=n)
n = midNodes.getByBoundingBox(yMin=width-tol, zMax=tol)
p.Set('BC-right-upper', nodes=n)
n = midNodes.getByBoundingBox(yMax=tol, zMin=thickness-tol)
p.Set('BC-left-lower', nodes=n)
n = midNodes.getByBoundingBox(yMax=tol, zMax=tol)
p.Set('BC-right-lower', nodes=n)

# Assign element type
p.Set('all', elements=p.elements)
p.setElementType(elemTypes=(mesh.ElemType(elemCode=C3D8R, elemLibrary=EXPLICIT, secondOrderAccuracy=OFF, hourglassControl=ENHANCED), ), regions=(p.sets['elastic']))
p.setElementType(elemTypes=(mesh.ElemType(elemCode=C3D8R, elemLibrary=EXPLICIT, secondOrderAccuracy=OFF, hourglassControl=ENHANCED), ), regions=(p.sets['damageable']))

# Assembly
instanceName = 'orphanMeshPart-1'
m.rootAssembly.DatumCsysByDefault(CARTESIAN)
m.rootAssembly.Instance(dependent=ON, name=instanceName, part=p)

# Step
m.ExplicitDynamicsStep(name='Step-1', previous='Initial', timePeriod=stepDuration, massScaling=((SEMI_AUTOMATIC, MODEL, THROUGHOUT_STEP, 0.0, 6.66667e-08, BELOW_MIN, 1, 0, 0.0, 0.0, 0, None), ), nlgeom=ON)

# Amplitude
m.SmoothStepAmplitude(name='Amp-1', timeSpan=STEP, data=((0.0, 0.0), (stepDuration, 1.0)))

# Boundary conditions
m.DisplacementBC(name='BC-fixz', createStepName='Step-1', region=m.rootAssembly.instances[instanceName].sets['allNodes'], u1=UNSET, u2=UNSET, u3=0.0, ur1=UNSET, ur2=UNSET, ur3=UNSET,
    amplitude='Amp-1', fixed=OFF, distributionType=UNIFORM, fieldName='', localCsys=None)
m.DisplacementBC(name='BC-pin', createStepName='Step-1', region=m.rootAssembly.instances[instanceName].sets['BC-pin'], u1=0.0, u2=0.0, u3=0.0, ur1=UNSET, ur2=UNSET, ur3=UNSET,
    amplitude='Amp-1', fixed=OFF, distributionType=UNIFORM, fieldName='', localCsys=None)
m.DisplacementBC(name='BC-clamp-left', createStepName='Step-1', region=m.rootAssembly.instances[instanceName].sets['BC-clamp-left'], u1=0.0, u2=0.0, u3=0.0, ur1=UNSET, ur2=UNSET, ur3=UNSET,
    amplitude='Amp-1', fixed=OFF, distributionType=UNIFORM, fieldName='', localCsys=None)
m.DisplacementBC(name='BC-clamp-right', createStepName='Step-1', region=m.rootAssembly.instances[instanceName].sets['BC-clamp-right'], u1=UNSET, u2=UNSET, u3=0.0, ur1=UNSET, ur2=UNSET, ur3=UNSET,
    amplitude='Amp-1', fixed=OFF, distributionType=UNIFORM, fieldName='', localCsys=None)
m.DisplacementBC(name='BC-loadApp', createStepName='Step-1', region=m.rootAssembly.instances[instanceName].sets['BC-loadApp'], u1=-0.01, u2=UNSET, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET,
    amplitude='Amp-1', fixed=OFF, distributionType=UNIFORM, fieldName='', localCsys=None)
m.Equation(name='LoadApp-1', terms=((-1.0, instanceName + '.BC-loadAppFollowers', 1), (1.0, instanceName + '.BC-loadApp', 1)))
m.Equation(name='y1', terms=((-1.0, instanceName + '.BC-left-upper', 1), (1.0, instanceName + '.BC-left-lower', 1)))
m.Equation(name='y2', terms=((-1.0, instanceName + '.BC-left-upper', 1), (1.0, instanceName + '.BC-right-lower', 1)))
m.Equation(name='y3', terms=((-1.0, instanceName + '.BC-left-upper', 1), (1.0, instanceName + '.BC-right-upper', 1)))

# Field output
m.FieldOutputRequest(createStepName='Step-1', numIntervals=50, name='damageable', region=m.rootAssembly.instances[instanceName].sets['all'], sectionPoints=DEFAULT, variables=('U', 'S', 'LE', 'SDV'))
m.fieldOutputRequests['F-Output-1'].setValues(numIntervals=50)

# History output - load app
m.HistoryOutputRequest(createStepName='Step-1', frequency=1, name='LoadApp', region=m.rootAssembly.instances[instanceName].sets['BC-loadApp'], sectionPoints=DEFAULT, variables=('U1', 'RF1'))
m.HistoryOutputRequest(createStepName='Step-1', frequency=1, name='plas', region=m.rootAssembly.instances[instanceName].sets['damageable'], sectionPoints=DEFAULT, variables=('SDV', ))

# Unsupported keywords
# *Depvar
m.keywordBlock.synchVersions(storeNodesAndElements=False)
pattern = re.compile(r'\*User Material,')
matches = [x for x in m.keywordBlock.sieBlocks if pattern.match(x)]
for match in matches:
    position = m.keywordBlock.sieBlocks.index(match)
    m.keywordBlock.insert(position,
        '*Depvar\n' +
        ' 26,\n' +
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
        ' 20, CDM_Plas13\n' +
        ' 21, CDM_Inel13\n' +
        ' 22, CDM_phi\n' +
        ' 23, CDM_gamma\n' +
        ' 24, CDM_Fm1\n' +
        ' 25, CDM_Fm2\n' +
        ' 26, CDM_Fm3\n' +
        '*Damping, alpha=1.d+4\n'+
        '*Characteristic Length, definition=USER, components=3\n'+
        '*Initial Conditions, Type=Solution\n'+
        ' orphanMeshPart-1.all,  0.d0,  0.d0,  0.d0,  0.d0,  0.d0,  0.d0,  0.d0,\n'+
        ' 0.d0,  0.d0,  -999,     1,  0.d0,  0.d0,  0.d0,  0.d0,\n'+
        ' 0.d0,  0.d0,  0.d0,  0.d0,  0.d0,  0.d0,  0.d0,  0.d0,\n'+
        ' 0.d0,  0.d0,  0.d0')

# pattern = re.compile(r'\*Node Output,')
# matches = [x for x in m.keywordBlock.sieBlocks if pattern.match(x)]
# if len(matches) == 1:
#     position = m.keywordBlock.sieBlocks.index(matches[0])
#     m.keywordBlock.insert(position,
#         '*Element Output, elset=orphanMeshPart-1.all\n' +
#         'SDV12, SDV20, SDV25 ')


# Create job
mdb.Job(model=modelName, name=jobName)
mdb.jobs[jobName].writeInput()
