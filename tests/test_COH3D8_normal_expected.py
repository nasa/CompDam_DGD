from utilities import get_enerElas

YT = 62.3
GYT = 0.277
length = 0.2

# Calculate elastic strain energy at area damage variable = 0.5, which is approx max value
enerElas_dmg_0p5 = get_enerElas(YT, GYT, length**2)

enerFrac = GYT * length**2
finalDisp = 2. * GYT / YT
finalDisp_tol = 0.005

parameters = {
	"results": [
		{
            "type": "max",
            "identifier":
                {
                    "symbol": "S33",
                    "elset": "COHESIVE",
                    "position": "Element 1 Int Point 1"
                },
            "referenceValue": YT,
            "tolerance": YT * 0.001  # 0.1% error
        },
        {
            "type": "disp_at_zero_y",
            "step": "Step-1",
            "identifier": [
                { # x
                    "symbol": "U3",
                    "nset": "Z+",
                    "position": "Node 5"
                },
                { # y
                    "symbol": "S33",
                    "elset": "COHESIVE",
                    "position": "Element 1 Int Point 1"
                }
            ],
            "window": [finalDisp*0.8, finalDisp*1.2],
            "zeroTol": 2e-4,  # Defines how close to zero the y value needs to be
            "referenceValue": finalDisp,  # Cohesive displacement-jump at complete failure
            "tolerance": finalDisp * finalDisp_tol
        },
        {
            "type": "max",
            "identifier":
                {
                    "symbol": "SDV_COH_dmg",
                    "elset": "COHESIVE",
                    "position": "Element 1 Int Point 1"
                },
            "referenceValue": 1.0,
            "tolerance": 0.0
        },
        {
            "type": "continuous",
            "identifier":
                {
                    "symbol": "S33",
                    "elset": "COHESIVE",
                    "position": "Element 1 Int Point 1"
                },
            "referenceValue": 0.0,
            "tolerance": 0.1
        },
        {
            "type": "max",
            "identifier": "Plastic dissipation: ALLPD for Whole Model",
            "referenceValue": enerFrac,  # Unrecoverable energy dissipation from fracture * fracture area: GYT*LC1*LC3
            "tolerance": enerFrac * 0.001  # 0.1% error
        },
        {
            "type": "finalValue",
            "identifier": "Plastic dissipation: ALLPD for Whole Model",
            "referenceValue": enerFrac,  # Unrecoverable energy dissipation from fracture * fracture area: GYT*LC1*LC3
            "tolerance": enerFrac * 0.001  # 0.1% error
        },
        {
            "type": "max",
            "identifier": "Strain energy: ALLSE for Whole Model",  # Recoverable strain energy
            "referenceValue": enerElas_dmg_0p5,  # Elastic strain energy * volume
            "tolerance": enerElas_dmg_0p5 * 0.001
        }
	]
}
