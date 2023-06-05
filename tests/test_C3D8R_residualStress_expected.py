sigma11_residual = -50.
sigma22_residual =  50.
sigma33_residual =   0.
tolerance = 0.00002

parameters = {
	"results": [
        {
            "type": "finalValue",
            "identifier":
                {
                    "symbol": "S11",
                    "elset": "ALL_ELEMS",
                    "position": "Element 1 Int Point 1"
                },
            "referenceValue": sigma11_residual,
            "tolerance": abs(sigma11_residual)*tolerance
        },
        {
            "type": "finalValue",
            "identifier":
                {
                    "symbol": "S22",
                    "elset": "ALL_ELEMS",
                    "position": "Element 1 Int Point 1"
                },
            "referenceValue": sigma22_residual,
            "tolerance": abs(sigma22_residual)*tolerance
        },
        {
            "type": "finalValue",
            "identifier":
                {
                    "symbol": "S33",
                    "elset": "ALL_ELEMS",
                    "position": "Element 1 Int Point 1"
                },
            "referenceValue": sigma33_residual,
            "tolerance": abs(sigma22_residual)*tolerance
        },
	]
}
