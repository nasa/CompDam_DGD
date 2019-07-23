aPL = 4.412E-10
nPL = 5.934
G12 = 5290.
enerPlas = aPL/G12*90.65**(nPL + 1.)*nPL/(nPL + 1.)*0.2**3
enerFrac = 0.788*0.2*0.2
enerTotal = enerPlas+enerFrac

parameters = {
"results": [
        {
            "type": "max",
            "step": "Step-7",
            "identifier":
                {
                    "symbol": "S12",
                    "elset": "ALL_ELEMS",
                    "position": "Element 1 Int Point 1"
                },
            "referenceValue": 90.65,
            "tolerance": 0.01
        },
        {
            "type": "max",
            "step": "Step-9",
            "identifier":
                {
                    "symbol": "S12",
                    "elset": "ALL_ELEMS",
                    "position": "Element 1 Int Point 1"
                },
            "referenceValue": 91.08,
            "tolerance": 0.01
        },
        {
            "type": "max",
            "step": "Step-11",
            "identifier":
                {
                    "symbol": "S12",
                    "elset": "ALL_ELEMS",
                    "position": "Element 1 Int Point 1"
                },
            "referenceValue": 79.96,
            "tolerance": 0.02
        },
        {
            "type": "max",
            "step": "Step-17",
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
            "step": "Step-17",
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
        },
        {
            "type": "tabular",
            "identifier": "Plastic dissipation: ALLPD for Whole Model",
            "referenceValue": [
                            (0.0, 0.0),
                            (0.10, aPL/G12*69.2368**(nPL + 1.)*nPL/(nPL + 1.)*0.2**3),  # End of step 1 (all plasticity)
                            (0.25, aPL/G12*79.045**(nPL + 1.)*nPL/(nPL + 1.)*0.2**3),  # End of step 3 (all plasticity)
                            (0.40, aPL/G12*85.5842**(nPL + 1.)*nPL/(nPL + 1.)*0.2**3),  # End of step 5 (all plasticity)
                            (0.55, aPL/G12*90.6525**(nPL + 1.)*nPL/(nPL + 1.)*0.2**3),  # End of step 7 (all plasticity)
                            ],
            "tolerance_percentage": 0.05
        },
        {
            "type": "finalValue",
            "identifier": "Plastic dissipation: ALLPD for Whole Model",
            "referenceValue": enerTotal,  # Total energy dissipation (plasticity + fracture)
            "tolerance": 0.02*enerTotal
        }
	]
}
