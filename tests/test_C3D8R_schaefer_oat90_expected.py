reference_values = [(0.000137549, 0.0312886),
                    (0.000629, 0.143),
                    (0.00103, 0.220102),
                    (0.0014, 0.26636)]
tolerance = [(a * 0.1, b * 0.1) for (a, b) in reference_values]
parameters = {
    "results": [
        {
            "type": "tabular",
            "identifier": [
                { # x
                    "symbol": "U2",
                    "nset": "INQUIRY"
                },
                { # y
                    "symbol": "RF2",
                    "nset": "INQUIRY"
                }
            ],
            "referenceValue": reference_values,
            "tolerance": tolerance
        }
    ]
}
