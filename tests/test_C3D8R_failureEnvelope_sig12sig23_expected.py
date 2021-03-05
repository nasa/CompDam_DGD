GSL = 0.788
length = 0.1

parameters = {
    "results": [
        {
            "type": "log_stress_at_failure_init",
            "step": "Step-1",
            "failureIndices": [
                {
                    "symbol": "SDV_CDM_FIM",
                    "elset": "ALL",
                    "position": "Element 1 Int Point 1"
                }
            ],
            "stressComponents": [
                {
                    "symbol": "S12",
                    "elset": "ALL",
                    "position": "Element 1 Int Point 1"
                },
                {
                    "symbol": "S23",
                    "elset": "ALL",
                    "position": "Element 1 Int Point 1"
                }
            ],
            "additionalIdentifiersToStore": [
                {
                    "symbol": "SDV_CDM_ALPHA",
                    "elset": "ALL",
                    "position": "Element 1 Int Point 1"
                }
            ]
        },
        {
            "type": "max",
            "step": "Step-1",
            "identifier": {
                "symbol": "SDV_CDM_FIm",
                "elset": "ALL",
                "position": "Element 1 Int Point 1"
            },
            "referenceValue": 1.0,
            "tolerance": 0.0
        },
        {
            "type": "max",
            "step": "Step-1",
            "identifier": {
                "symbol": "SDV_CDM_d1T",
                "elset": "ALL",
                "position": "Element 1 Int Point 1"
            },
            "referenceValue": 0.0,
            "tolerance": 0.0
        },
        {
            "type": "max",
            "step": "Step-1",
            "identifier": {
                "symbol": "SDV_CDM_d1C",
                "elset": "ALL",
                "position": "Element 1 Int Point 1"
            },
            "referenceValue": 0.0,
            "tolerance": 0.0
        },
        {
            "type": "max",
            "identifier": "Plastic dissipation: ALLPD for Whole Model",
            "referenceValue": GSL * (length * length),  # Unrecoverable energy dissipation from fracture * fracture area: GSL*LC1*LC3
            "tolerance": GSL * (length * length) * 0.01
        }
    ]
}
