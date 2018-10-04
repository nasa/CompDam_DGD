reference_values = [(0.0000123350125495656, 0.108519895263409),
                    (0.002315396, 20.3694371436868),
                    (0.00473865803221036, 41.6163221302831),
                    (0.00711880771429827, 61.7194939233583),
                    (0.00971695556460476, 80.5043944400537),
                    (0.0109375645469152, 87.6431991939856),
                    ]
length = 0.1
area = length ** 2
tolerance = [(a * 0.05, b * 0.05) for (a, b) in reference_values]
parameters = {
	"results": [
        {
            "type": "tabular",
            "step": "Step-1",
            "identifier": [
                {   "av_id": "n3_u",
                    "symbol": "U2",
                    "nset": "Y+",
                    "position": "Node 3"
                },
                {   "av_id": "n3_rf",
                    "symbol": "RF2",
                    "nset": "Y+",
                    "position": "Node 3"
                },
                {   "av_id": "n4_rf",
                    "symbol": "RF2",
                    "nset": "Y+",
                    "position": "Node 4"
                }
            ],
            "x_eval_statement": "d['n3_u'] / {}".format(length),
            "y_eval_statement": "(d['n3_rf'] + d['n4_rf']) / {}".format(area),          
            "referenceValue": reference_values,
            "tolerance": tolerance
        }
	]
}