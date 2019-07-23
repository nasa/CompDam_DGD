E1 = 171420.
strain = 0.02
enerElas = 0.5*strain**2*E1*0.5*0.5*0.5

parameters = {
    "results": [
        {
            "type": "slope",
            "step": "Step-1",
            "identifier": [
                { # x
                    "symbol": "LE11",
                    "elset": "ALL",
                    "position": "Element 1 Int Point 1"
                },
                { # y
                    "symbol": "S11",
                    "elset": "ALL",
                    "position": "Element 1 Int Point 1"
                }
            ],
            "window": [0.0001, 0.005],  # [min, max] in x
            "referenceValue": E1,
            "tolerance": 0.01*E1 # 1%
        },
        {
            "type": "max",
            "identifier":
                {
                    "symbol": "SDV_CDM_d2",
                    "elset": "ALL",
                    "position": "Element 1 Int Point 1"
                },
            "referenceValue": 0.0
        },
        {
            "type": "max",
            "identifier":
                {
                    "symbol": "SDV_CDM_d1T",
                    "elset": "ALL",
                    "position": "Element 1 Int Point 1"
                },
            "referenceValue": 0.0
        },
        {
            "type": "max",
            "identifier":
                {
                    "symbol": "SDV_CDM_d1C",
                    "elset": "ALL",
                    "position": "Element 1 Int Point 1"
                },
            "referenceValue": 0.0
        },
        {
            "type": "max",
            "identifier":
                {
                    "symbol": "SDV_CDM_FIm",
                    "elset": "ALL",
                    "position": "Element 1 Int Point 1"
                },
            "referenceValue": 0.0
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
            "tolerance": 0.1
        },
        {
            "type": "max",
            "identifier": "Strain energy: ALLSE for Whole Model",   # Recoverable strain energy
            "referenceValue": enerElas,   # Elastic strain energy * volume
            "tolerance": 0.03*enerElas
        }
    ]
}
