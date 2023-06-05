aPL = 4.412E-10
nPL = 5.934
G12 = 5290.
GSL = 0.788
length = 0.2

plastic_energy_dissipated = lambda stress : aPL / G12 * stress**(nPL + 1.) * nPL / (nPL + 1.) * length**3

enerPlas = plastic_energy_dissipated(90.65)
enerFrac = GSL * (length * length)
enerTotal = enerPlas + enerFrac

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
            "tolerance": 0.04
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
                            (0.10, plastic_energy_dissipated(69.2368)),  # End of step 1 (all plasticity)
                            (0.25, plastic_energy_dissipated(79.0450)),  # End of step 3 (all plasticity)
                            (0.40, plastic_energy_dissipated(85.5842)),  # End of step 5 (all plasticity)
                            (0.55, plastic_energy_dissipated(90.6525)),  # End of step 7 (all plasticity)
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
