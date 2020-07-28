staticLoad = 2.352

parameters = {
    "results": [
		{
            "type": "max",
			"step": "Load",
            "identifier": {
                "symbol": "RF3",
                "nset": "T_NODE",
                "position": "Node 9999998"
            },
            "referenceValue": staticLoad,
            "tolerance": staticLoad * 0.01
        },
        {
            "type": "max",
            "identifier": "Plastic dissipation: ALLPD in ELSET COH",
            "referenceValue": 1.605,  # Unrecoverable energy dissipation from fracture * fracture area: GIc_init * delta_a * width
            "tolerance": 0.05
        },
    ]
}
