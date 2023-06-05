import math as m

XT = 2326.2   # Fiber tensile strength
E1 = 171420.  # Young's modulus
nu12 = 0.32   # Poisson ratio
GXT = 133.3   # Fiber tensile fracture toughness
fXT = 0.75    # Fiber tensile strength ratio
fGXT = 0.375  # Fiber tensile fracture toughness ratio
Lc = 0.5      # Fiber-direction element size

# Green-Lagrange strains of interest
epsilon_initiation = XT/E1
epsilon_final_A = 2.*GXT*fGXT / (XT*fXT*Lc)  # 'A' represents the first superposed softening law, scaled by n and m
epsilon_final_B = 2.*GXT*(1.-fGXT) / (XT*(1.-fXT)*Lc)  # 'B' represents the second superposed softening law
epsilon_inflection = max(epsilon_initiation, min(epsilon_final_A, epsilon_final_B))
epsilon_final      = max(epsilon_final_A, epsilon_final_B)

if epsilon_final_B > epsilon_final_A:
    sigma_inflection = XT*(1. - fXT)*(epsilon_final - epsilon_inflection)/(epsilon_final - epsilon_initiation)
else:
    sigma_inflection = XT*fXT*(epsilon_final - epsilon_inflection)/(epsilon_final - epsilon_initiation)
    
# Displacements of interest
green2stretch = lambda green_strain : m.sqrt(2. * green_strain + 1.)
displacement_inflection = (green2stretch(epsilon_inflection)-1)*Lc
displacement_final = (green2stretch(epsilon_final)-1)*Lc

# Cauchy stresses
Cauchy_initation = E1 * epsilon_initiation / ((green2stretch(-nu12*epsilon_initiation)**2) / green2stretch(epsilon_initiation))


parameters = {
    "results": [
        {
            "type": "max",
            "identifier":
                {
                    "symbol": "S11",
                    "elset": "ALL",
                    "position": "Element 1 Int Point 1"
                },
            "referenceValue": Cauchy_initation,
            "tolerance": Cauchy_initation * 0.001  # 0.1% error
        },
        {
            "type": "xy_infl_pt",
            "identifier": [
                { # x
                    "symbol": "U1",
                    "nset": "LOADAPP"
                },
                { # y
                    "symbol": "S11",
                    "elset": "ALL",
                    "position": "Element 1 Int Point 1"
                }
            ],
            "window": [displacement_inflection*0.9, displacement_inflection*1.1],
            "referenceValue": [displacement_inflection, sigma_inflection],
            "tolerance": [displacement_inflection*0.001, sigma_inflection*0.15]
        },
        {
            "type": "disp_at_zero_y",
            "identifier": [
                { # x
                    "symbol": "U1",
                    "nset": "LOADAPP"
                },
                { # y
                    "symbol": "S11",
                    "elset": "ALL",
                    "position": "Element 1 Int Point 1"
                }
            ],
            "window": [displacement_final*0.9, displacement_final*1.1],  # [min, max]
            "zeroTol": XT*0.001,  # Defines how close to zero the y value needs to be
            "referenceValue": displacement_final,
            "tolerance": displacement_final*0.01
        },
        {
            "type": "max",
            "identifier":
                {
                    "symbol": "SDV_CDM_d1T",
                    "elset": "ALL",
                    "position": "Element 1 Int Point 1"
                },
            "referenceValue": 1.0,
            "tolerance": 0.000001
        },
        {
            "type": "continuous",
            "identifier":
                {
                    "symbol": "S11",
                    "elset": "ALL",
                    "position": "Element 1 Int Point 1"
                },
            "referenceValue": XT/100000.,
            "tolerance": XT/1000.
        }
    ]
}
