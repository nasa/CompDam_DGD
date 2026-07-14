import os, sys

# Make trilinear-tsl utilities available for import
compdam_source_dir = os.environ.get('COMPDAM_SOURCE_DIR')
if compdam_source_dir is None:
    this_file_dir = os.path.dirname(os.path.abspath(__file__))
    utilities_dir = os.path.join(this_file_dir, '..', 'utilities', 'trilinear-cdm')
else:
    utilities_dir = os.path.join(compdam_source_dir, 'utilities', 'trilinear-cdm')
sys.path.append(utilities_dir)

# Import the trilinear law CDM implementation
from trilinearCDM import cauchy_stress_at_initiation, cauchy_stress_at_inflection, displacement_at_inflection, displacement_at_final

kw = dict(
    XC = -1200.1,   # Fiber compressive strength
    E1 = 171420.,  # Young's modulus
    E2 = 9080,
    G12 = 5290,
    nu12 = 0.32,   # Poisson ratio
    nu23 = 0.52,
    GXC = 61.0,   # Fiber compressive fracture toughness
    fXC = 0.614,    # Fiber compressive strength ratio
    fGXC = 0.432,  # Fiber compressive fracture toughness ratio
    Lc = 0.5      # Fiber-direction element size
)

cauchy_initiation = cauchy_stress_at_initiation(**kw)
displacement_inflection = displacement_at_inflection(**kw)
cauchy_inflection11 = cauchy_stress_at_inflection(**kw)
displacement_final = displacement_at_final(**kw)


parameters = {
	"results": [
        {
            "type": "min",
            "step": "Step-1",
            "identifier": {
                "symbol": "S11",
                "elset": "ALL",
                "position": "Element 1 Int Point 1"
            },
            "referenceValue": cauchy_initiation,
            "tolerance": -1*cauchy_initiation * 0.001  # 0.1% error
        },
        {
            "type": "xy_infl_pt",
            "step": "Step-1",
            "identifier": [
                { # x
                    "symbol": "U1",
                    "nset": "LOADAPP"
                },
                { # y
                    "symbol": "S11",
                    "elset": "ALL",
                    "position": "Element 1 Int Point 1"
                }
            ],
            "window": [displacement_inflection*1.1, displacement_inflection*0.9],
            "referenceValue": [displacement_inflection, cauchy_inflection11],
            "tolerance": [-1*displacement_inflection*0.001, -1*cauchy_inflection11*0.005]
        },
        {
            "type": "disp_at_zero_y",
            "step": "Step-1",
            "identifier": [
                { # x
                    "symbol": "U1",
                    "nset": "LOADAPP"
                },
                { # y
                    "symbol": "S11",
                    "elset": "ALL",
                    "position": "Element 1 Int Point 1"
                }
            ],
            "window": [displacement_final*1.5, displacement_final*0.9],  # [min, max]
            "zeroTol": -1*kw["XC"]*0.001,  # Defines how close to zero the y value needs to be
            "referenceValue": displacement_final,
            "tolerance": -1*displacement_final*0.01
        },
        {
            "type": "max",
            "step": "Step-1",
            "identifier": {
                "symbol": "SDV_CDM_D2",
                "elset": "ALL",
                "position": "Element 1 Int Point 1"
            },
            "referenceValue": 0.0,
            "tolerance": 0.0
        },
        {
            "type": "finalValue",
            "identifier": {
                "symbol": "SDV_CDM_STATUS",
                "elset": "ALL",
                "position": "Element 1 Int Point 1"
            },
            "referenceValue": 0.0,
            "tolerance": 0.0
        },
        {
            "type": "continuous",
            "identifier": {
                "symbol": "S11",
                "elset": "ALL",
                "position": "Element 1 Int Point 1"
            },
            "referenceValue": 0.,
            "tolerance": -1*kw["XC"]/1000.
        }
	]
}
