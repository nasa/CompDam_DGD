import math as m

E1 = 171420.
E2 = 9080.
nu12 = 0.32
nu23 = 0.52

length = 0.2
applied_engineering_strain = 0.02

nu21 = nu12 * E2 / E1

log2green = lambda log_strain : 0.5 * (m.exp(2. * log_strain) - 1.)
eng2green = lambda engineering_strain : 0.5 * ((1. + engineering_strain)**2 - 1.)
green2stretch = lambda green_strain : m.sqrt(2. * green_strain + 1.)

enerElas = 0.5 * eng2green(applied_engineering_strain)**2 * E2 * (length*length*length)

# Strain range for modulus check
log_strain_1 = 0.0001
log_strain_2 = 0.005

# 2nd Piola-Kirchhoff stresses for strain range
PK2_stress_2 = E2 * log2green(log_strain_2)
PK2_stress_1 = E2 * log2green(log_strain_1)

# Cauchy stresses for strain range
Cauchy_2 = PK2_stress_2 / ((green2stretch(-nu21*log2green(log_strain_2))*green2stretch(-nu23*log2green(log_strain_2))) / green2stretch(log2green(log_strain_2)))
Cauchy_1 = PK2_stress_1 / ((green2stretch(-nu21*log2green(log_strain_1))*green2stretch(-nu23*log2green(log_strain_1))) / green2stretch(log2green(log_strain_1)))

parameters = {
	"results": [
        {
            "type": "slope",
            "step": "Step-1",
            "identifier": [
                { # x, Logarithmic strain
                    "symbol": "LE22",
                    "elset": "ALL_ELEMS",
                    "position": "Element 1 Int Point 1"
                },
                { # y, Cauchy stress
                    "symbol": "S22",
                    "elset": "ALL_ELEMS",
                    "position": "Element 1 Int Point 1"
                }
            ],
            "window": [log_strain_1, log_strain_2],  # [min, max] in x
            "referenceValue": (Cauchy_2 - Cauchy_1) / (log_strain_2 - log_strain_1),
            "tolerance": (Cauchy_2 - Cauchy_1) / (log_strain_2 - log_strain_1) * 0.005  # 0.5% tolerance
        },
        {
            "type": "max",
            "identifier":
                {
                    "symbol": "SDV_CDM_d2",
                    "elset": "ALL_ELEMS",
                    "position": "Element 1 Int Point 1"
                },
            "referenceValue": 0.0
        },
        {
            "type": "max",
            "identifier":
                {
                    "symbol": "SDV_CDM_d1T",
                    "elset": "ALL_ELEMS",
                    "position": "Element 1 Int Point 1"
                },
            "referenceValue": 0.0
        },
        {
            "type": "max",
            "identifier":
                {
                    "symbol": "SDV_CDM_d1C",
                    "elset": "ALL_ELEMS",
                    "position": "Element 1 Int Point 1"
                },
            "referenceValue": 0.0
        },
        {
            "type": "max",
            "identifier":
                {
                    "symbol": "SDV_CDM_FIm",
                    "elset": "ALL_ELEMS",
                    "position": "Element 1 Int Point 1"
                },
            "referenceValue": 0.0
        },
        {
            "type": "continuous",
            "identifier":
                {
                    "symbol": "S22",
                    "elset": "ALL_ELEMS",
                    "position": "Element 1 Int Point 1"
                },
            "referenceValue": 0.0,
            "tolerance": 0.1
        },
        {
            "type": "max",
            "identifier": "Strain energy: ALLSE for Whole Model",  # Recoverable strain energy
            "referenceValue": enerElas,  # Elastic strain energy * volume
            "tolerance": enerElas*0.01  # 1% tolerance
        }
	]
}
