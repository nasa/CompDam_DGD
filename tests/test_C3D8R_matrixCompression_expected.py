YC = 199.8  # matrix compression strength
length = 0.15
thickness = 0.1

parameters = {
    "results": [
        {
            "type": "min",
            "identifier":
                {
                    "symbol": "S22",
                    "elset": "ALL",
                    "position": "Element 1 Int Point 1"
                },
            "referenceValue": -YC,
            "tolerance": YC * 0.005
        },
        {
            "type": "max",
            "identifier":
                {
                    "symbol": "SDV_CDM_d2",
                    "elset": "ALL",
                    "position": "Element 1 Int Point 1"
                },
            "referenceValue": 1.0,
            "tolerance": 0.0
        },
        {
            "type": "max",
            "identifier":
                {
                    "symbol": "SDV_CDM_d1T",
                    "elset": "ALL",
                    "position": "Element 1 Int Point 1"
                },
            "referenceValue": 0.0,
            "tolerance": 0.0
        },
        {
            "type": "max",
            "identifier":
                {
                    "symbol": "SDV_CDM_d1C",
                    "elset": "ALL",
                    "position": "Element 1 Int Point 1"
                },
            "referenceValue": 0.0,
            "tolerance": 0.0
        },
        {
            "type": "finalValue",
            "identifier":
                {
                    "symbol": "SDV_CDM_alpha",
                    "elset": "ALL",
                    "position": "Element 1 Int Point 1"
                },
            "referenceValue": 53.0,
            "tolerance": 0.4
        },
        {
            "type": "max",
            "identifier": "Plastic dissipation: ALLPD for Whole Model",
            "referenceValue": 0.788 * (length * thickness),  # Unrecoverable energy dissipation from fracture * fracture area: GSL*LC1*LC3
            "tolerance": 0.788 * (length * thickness) * 0.01
        },
        {
            "type": "continuous",
            "identifier":
                {
                    "symbol": "S22",
                    "elset": "ALL",
                    "position": "Element 1 Int Point 1"
                },
            "referenceValue": 0.0,
            "tolerance": 2.0
        }
    ]
}
