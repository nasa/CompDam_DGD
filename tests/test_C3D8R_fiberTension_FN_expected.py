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
    XT = 2326.2,   # Fiber tensile strength
    E1 = 171420.,  # Young's modulus
    E2 = 9080,
    G12 = 5290,
    nu12 = 0.32,   # Poisson ratio
    nu23 = 0.52,
    GXT = 205.0,   # Fiber tensile fracture toughness
    fXT = 0.63,    # Fiber tensile strength ratio
    fGXT = 0.469,  # Fiber tensile fracture toughness ratio
    Lc = 0.5      # Fiber-direction element size
)

cauchy_initiation = cauchy_stress_at_initiation(**kw)
displacement_inflection = displacement_at_inflection(**kw)
cauchy_inflection11 = cauchy_stress_at_inflection(**kw)
displacement_final = displacement_at_final(**kw)


parameters = {
    "results": [
        {
            "type": "max",
            "identifier":
                {
                    "symbol": "S11",
                    "elset": "ALL",
                    "position": "Element 1 Int Point 1"
                },
            "referenceValue": cauchy_initiation,
            "tolerance": cauchy_initiation * 0.005  # 0.5% error
        },
        {
            "type": "xy_infl_pt",
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
            "window": [displacement_inflection*0.9, displacement_inflection*1.1],
            "referenceValue": [displacement_inflection, cauchy_inflection11],
            "tolerance": [displacement_inflection*0.001, cauchy_inflection11*0.005]
        },
        {
            "type": "disp_at_zero_y",
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
            "window": [displacement_final*0.7, displacement_final*1.1],  # [min, max]
            "zeroTol": kw["XT"]*0.001,  # Defines how close to zero the y value needs to be
            "referenceValue": displacement_final,
            "tolerance": displacement_final*0.01
        },
        {
            "type": "max",
            "identifier":
                {
                    "symbol": "SDV_CDM_d1T",
                    "elset": "ALL",
                    "position": "Element 1 Int Point 1"
                },
            "referenceValue": 1.0,
            "tolerance": 0.000001
        },
        {
            "type": "continuous",
            "identifier":
                {
                    "symbol": "S11",
                    "elset": "ALL",
                    "position": "Element 1 Int Point 1"
                },
            "referenceValue": kw["XT"]/100000.,
            "tolerance": kw["XT"]/1000.
        }
    ]
}
