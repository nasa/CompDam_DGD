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
            "tolerance": staticLoad * 0.005
        },
        {
            "type": "max",
            "identifier": "Plastic dissipation: ALLPD for Whole Model",
            "referenceValue": 1.63141,  # Unrecoverable energy dissipation from fracture * fracture area: GYT * delta_a * width
            "tolerance": 0.01
        },
    ]
}
