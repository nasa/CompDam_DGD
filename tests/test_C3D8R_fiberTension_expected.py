import math as m

XT = 2326.2  # Fiber tensile strength
E1 = 171420.  # Young's modulus
nu12 = 0.32

failure_strain = XT / E1  # in terms of Green-Lagrange strain

green2stretch = lambda green_strain : m.sqrt(2. * green_strain + 1.)

# Cauchy stresses for strain range
Cauchy_at_initation = E1 * failure_strain / ((green2stretch(-nu12*failure_strain)**2) / green2stretch(failure_strain))

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
            "referenceValue": Cauchy_at_initation,
            "tolerance": Cauchy_at_initation * 0.001  # 0.1% error
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
            "window": [0.05, 0.06],
            "referenceValue": [0.05730, 465.24],
            "tolerance": [0.003, 65]
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
            "window": [0.15, 0.4],  # [min, max]
            "zeroTol": 0.00623,  # Defines how close to zero the y value needs to be
            "referenceValue": 0.2865,
            "tolerance": 0.055
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
            "type": "max",
            "identifier":
                {
                    "symbol": "SDV_CDM_d2",
                    "elset": "ALL",
                    "position": "Element 1 Int Point 1"
                },
            "referenceValue": 0.0,
            "tolerance": 0.0
        },
        {
            "type": "continuous",
            "identifier":
                {
                    "symbol": "S11",
                    "elset": "ALL",
                    "position": "Element 1 Int Point 1"
                },
            "referenceValue": 0.0,
            "tolerance": 0.2
        }
    ]
}
