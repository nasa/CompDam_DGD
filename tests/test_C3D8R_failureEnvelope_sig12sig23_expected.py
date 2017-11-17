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
        }
    ]
}
