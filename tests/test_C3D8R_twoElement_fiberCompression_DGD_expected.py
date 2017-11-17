# Variable results
crushStress=-9.55

parameters = {
    "results": [
        {
            "type": "min",
            "step": "Step-1",
            "identifier": {
                    "symbol": "RF1",
                    "nset": "BC-LOADAPP"
                },
            "referenceValue": -36.003,
            "tolerance": 0.72  # 2%
        },
        {
            "type": "finalValue",
            "step": "Step-1",
            "identifier": {
                "symbol": "RF1",
                "nset": "BC-LOADAPP",
            },
            "referenceValue": crushStress,
            "tolerance": 0.02*abs(crushStress)  # 2%
        }
    ]
}