import math

applied_compression = -40.0
coefficient_of_friction = 0.3

SL = 92.3  # longitudinal shear strength
YC = 199.8  # matrix compressive strength
G2C = 0.788  # Mode II matrix fracture toughness
alpha0_DGD = 0.9080748  # matrix crack orientation due to pure matrix compression failure
length = 0.1  # element edge length

eta_L = -SL * math.cos(2*alpha0_DGD) / (YC * math.cos(alpha0_DGD) * math.cos(alpha0_DGD))

parameters = {
"results": [
        {

            "type": "min",
            "step": "Step-1",
            "identifier":
                {
                    "symbol": "S22",
                    "elset": "ALL_ELEMS",
                    "position": "Element 1 Int Point 1"
                },
            "referenceValue": applied_compression,
            "tolerance": abs(applied_compression) * 0.005
        },
        {
            "type": "max",
            "step": "Step-2",
            "identifier":
                {
                    "symbol": "S12",
                    "elset": "ALL_ELEMS",
                    "position": "Element 1 Int Point 1"
                },
            "referenceValue": SL - eta_L * applied_compression,
            "tolerance": (SL - eta_L * applied_compression) * 0.005
        },
        {
            "type": "max",
            "step": "Step-3",
            "identifier":
                {
                    "symbol": "S12",
                    "elset": "ALL_ELEMS",
                    "position": "Element 1 Int Point 1"
                },
            "referenceValue": -applied_compression * coefficient_of_friction,
            "tolerance": abs(applied_compression) * coefficient_of_friction * 0.005
        },
        {
            "type": "min",
            "step": "Step-3",
            "identifier":
                {
                    "symbol": "S12",
                    "elset": "ALL_ELEMS",
                    "position": "Element 1 Int Point 1"
                },
            "referenceValue": applied_compression * coefficient_of_friction,
            "tolerance": abs(applied_compression) * coefficient_of_friction * 0.005
        },
        {
            "type": "max",
            "identifier": "Plastic dissipation: ALLPD for Whole Model",
            "referenceValue": G2C * (length * length),  # Unrecoverable energy dissipation from fracture * fracture area: GSL*LC1*LC3
            "tolerance": G2C * (length * length) * 0.005
        },
        {
            "type": "max",
            "step": "Step-3",
            "identifier":
                {
                    "symbol": "SDV_CDM_d1T",
                    "elset": "ALL_ELEMS",
                    "position": "Element 1 Int Point 1"
                },
            "referenceValue": 0.0,
            "tolerance": 0.0
        },
        {
            "type": "max",
            "step": "Step-3",
            "identifier":
                {
                    "symbol": "SDV_CDM_d1C",
                    "elset": "ALL_ELEMS",
                    "position": "Element 1 Int Point 1"
                },
            "referenceValue": 0.0,
            "tolerance": 0.0
        },
        {
            "type": "continuous",
            "identifier":
                {
                    "symbol": "S12",
                    "elset": "ALL_ELEMS",
                    "position": "Element 1 Int Point 1"
                },
            "referenceValue": 0.0,
            "tolerance": 0.1
        }
	]
}
