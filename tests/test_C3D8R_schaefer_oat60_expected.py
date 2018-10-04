reference_values = [(0.000000,      0.000000),
                    (0.000217,      2.245748),
                    (0.001460,     14.883391),
                    (0.004090,     41.966763),
                    (0.007958,     77.524742),
                    (0.012524,    103.247925),
                    (0.017066,    116.964104),
                    (0.020941,    124.527130),
                    (0.023559,    129.381104),
                    (0.024789,    131.595917),
                    (0.025000,    132.018784), ]
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