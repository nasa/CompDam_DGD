# (strain, stress) value pairs 
reference_values = [(0.000000,      0.000000),
                    (0.000217,      2.031003),
                    (0.001449,     13.514810),
                    (0.004098,     38.584564),
                    (0.007939,     71.311012),
                    (0.012510,     99.755920),
                    (0.017078,    115.675209),
                    (0.020937,    124.543907),
                    (0.023561,    129.376526),
                    (0.024788,    131.681854),
                    (0.025000,    132.192749),]
length = 0.1
area = length * length
tolerance = [(a * 0.05, b * 0.05) for (a, b) in reference_values]
parameters = {
    "results": [
        {
            "type": "tabular",
            "identifier": [
                { 
                    "label": "u_n5",
                    "symbol": "U1",
                    "position": "Node 5",
                    "nset": "RIGHTEDGE"
                },
                { 
                    "label": "rf_n5",
                    "symbol": "RF1",
                    "position": "Node 5",
                    "nset": "RIGHTEDGE"
                },
                { 
                    "label": "rf_n6",
                    "symbol": "RF1",
                    "position": "Node 6",
                    "nset": "RIGHTEDGE"
                },
                { 
                    "label": "rf_n7",
                    "symbol": "RF1",
                    "position": "Node 7",
                    "nset": "RIGHTEDGE"
                },
                { 
                    "label": "rf_n8",
                    "symbol": "RF1",
                    "position": "Node 8",
                    "nset": "RIGHTEDGE"
                },
            ],
            "xEvalStatement": "d['u_n5'] / {}".format(length),
            "yEvalStatement": "(d['rf_n5'] + d['rf_n6'] + d['rf_n7'] + d['rf_n8']) / {}".format(area),          
            "referenceValue": reference_values,
            "tolerance": tolerance
        }
    ]
}
