reference_values = [(0.000000,      0.000000),
                    (0.000058,      1.315872),
                    (0.000427,      9.549489),
                    (0.001313,     29.437089),
                    (0.002826,     63.485028),
                    (0.004987,    107.413523),
                    (0.007755,    142.915494),
                    (0.011031,    166.504341),
                    (0.014675,    185.256466),
                    (0.018510,    200.025433),
                    (0.024228,    217.336078), ]
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