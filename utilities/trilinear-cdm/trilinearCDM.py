import numpy as np
try:
    import matplotlib.pyplot as plt
except ImportError:
    plt = None

def generic_variable_names_from_dict(**kwargs):
    """
    Helper function so that inputs can be provided using either tension or compression naming
    conventions.
    """
    if 'strength' in kwargs:
        strength = kwargs['strength']
    elif 'XT' in kwargs:
        strength = kwargs['XT']
    elif 'XC' in kwargs:
        strength = kwargs['XC']
    else:
        strength = None

    if 'toughness' in kwargs:
        GX = kwargs['toughness']
    elif 'GX' in kwargs:
        GX = kwargs['GX']
    elif 'GXT' in kwargs:
        GX = kwargs['GXT']
    elif 'GXC' in kwargs:
        GX = kwargs['GXC']
    else:
        GX = None

    if 'n' in kwargs:
        n = kwargs['n']
    elif 'fXT' in kwargs:
        n = kwargs['fXT']
    elif 'fXC' in kwargs:
        n = kwargs['fXC']
    else:
        n = None

    if 'm' in kwargs:
        m = kwargs['m']
    elif 'fGXT' in kwargs:
        m = kwargs['fGXT']
    elif 'fGXC' in kwargs:
        m = kwargs['fGXC']
    else:
        m = None

    E1 = kwargs.get('E1', None)
    E2 = kwargs.get('E2', None)
    G12 = kwargs.get('G12', None)
    nu12 = kwargs.get('nu12', None)
    nu23 = kwargs.get('nu23', None)
    Lc = kwargs.get('Lc', None)

    return dict(
        strength = strength,
        GX = GX,
        n = n,
        m = m,
        E1 = E1,
        E2 = E2,
        G12 = G12,
        nu12 = nu12,
        nu23 = nu23,
        Lc = Lc
    )

def stretch_from_gl_strain(green_strain):
    return np.sqrt(2. * green_strain + 1.)

def gl_strain_from_stretch(stretch):
    return 0.5*(stretch**2 - 1.)

def stiffness6x6(E1, E2, G12, nu12, nu23, d1=0, d2=0, d3=0):
    E3 = E2
    nu13 = nu12
    G13 = G12
    G23 = E2/2/(1+nu23)
    nu21 = nu12*E2/E1
    nu31 = nu13*E3/E1
    nu32 = nu23*E3/E2

    delta = 1/(1 - nu12*nu21*(1-d1)*(1-d2) - nu23*nu32*(1-d2)*(1-d3) - nu31*nu13*(1-d1)*(1-d3) - 2*nu12*nu23*nu31*(1-d1)*(1-d2)*(1-d3))
    Stiff = np.zeros((6,6))

    Stiff[0,0] = delta*(1-d1)*(1 - nu23*nu32*(1-d2)*(1-d3))*E1
    Stiff[1,1] = delta*(1-d2)*(1 - nu13*nu31*(1-d1)*(1-d3))*E2
    Stiff[0,1] = delta*(1-d1)*(1-d2)*(nu12 + nu32*nu13*(1-d3))*E2
    Stiff[1,0] = Stiff[0,1]
    Stiff[2,2] = delta*(1-d3)*(1 - nu12*nu21*(1-d1)*(1-d2))*E3
    Stiff[3,3] = (1-d1)*(1-d2)*G12
    Stiff[0,2] = delta*(1-d1)*(1-d3)*(nu13+nu12*nu23*(1-d2))*E3
    Stiff[2,0] = Stiff[0,2]
    Stiff[1,2] = delta*(1-d2)*(1-d3)*(nu23+nu21*nu13*(1-d1))*E3
    Stiff[2,1] = Stiff[1,2]
    Stiff[4,4] = (1-d2)*(1-d3)*G23
    Stiff[5,5] = (1-d1)*(1-d3)*G13
    return Stiff

def gl_strain_cauchy_stress_from_stretch(stretch, stiff):
    if isinstance(stretch, (int, float)):
        stretch = np.array([stretch])
    eps11 = np.zeros_like(stretch)
    cauchy = np.zeros((3,3,len(stretch)))
    for i, s in enumerate(stretch):
        F = np.array([[s, 0., 0.], [0., 1., 0.], [0., 0., 1.]])
        eps = 0.5*(np.matmul(F.T, F) - np.eye(3))
        eps_vec = np.array([eps[0,0], eps[1,1], eps[2,2], 2*eps[0,1], 2*eps[1,2], 2*eps[0,2]])
        pk2_vec = np.matmul(stiff, eps_vec)
        pk2 = np.array([[pk2_vec[0], pk2_vec[3], pk2_vec[5]], [pk2_vec[3], pk2_vec[1], pk2_vec[4]], [pk2_vec[5], pk2_vec[4], pk2_vec[2]]])
        cauchy[:,:,i] = np.matmul(np.matmul(F.T, pk2), F) / np.linalg.det(F)
        eps11[i] = eps[0,0]
    return (eps11, cauchy)

def cauchy_stress_at_initiation(**kwargs):
    args = generic_variable_names_from_dict(**kwargs)
    strength = args['strength']
    E1 = args['E1']
    nu12 = args['nu12']
    if None in [strength, E1, nu12]:
        raise KeyError("cauchy_stress_at_initiation requires strength, E1, and nu12 to be defined in kwargs.")
    epsilon_initiation = strength/E1
    green2stretch = lambda green_strain : np.sqrt(2. * green_strain + 1.)
    return E1 * epsilon_initiation / ((green2stretch(-nu12*epsilon_initiation)**2) / green2stretch(epsilon_initiation))

def cauchy_stress_at_inflection(**kwargs):
    args = generic_variable_names_from_dict(**kwargs)
    strength = args['strength']
    E1 = args['E1']
    E2 = args['E2']
    G12 = args['G12']
    nu12 = args['nu12']
    nu23 = args['nu23']
    GX = args['GX']
    n = args['n']
    m = args['m']
    Lc = args['Lc']
    if None in [strength, E1, E2, G12, nu12, nu23, GX, n, m, Lc]:
        raise KeyError("cauchy_stress_at_inflection requires strength, E1, E2, G12, nu12, nu23, GX, n, m, and Lc to be defined in kwargs.")
    epsilon_initiation = strength/E1
    epsilon_final_A = 2.*GX*m / (strength*n*Lc)  # 'A' represents the first superposed softening law, scaled by n and m
    epsilon_final_B = 2.*GX*(1.-m) / (strength*(1.-n)*Lc)  # 'B' represents the second superposed softening law
    if strength > 0:  # Tension
        epsilon_inflection = max(epsilon_initiation, min(epsilon_final_A, epsilon_final_B))
        epsilon_final      = max(epsilon_final_A, epsilon_final_B)
    else:  # Compression
        epsilon_inflection = min(epsilon_initiation, max(epsilon_final_A, epsilon_final_B))
        epsilon_final      = min(epsilon_final_A, epsilon_final_B)

    # Stress at inflection point - input value (will differ from output when strains are large)
    if abs(epsilon_final_B) > abs(epsilon_final_A):
        sigma_inflection = strength*(1. - n)*(epsilon_final - epsilon_inflection)/(epsilon_final - epsilon_initiation)
    else:
        sigma_inflection = strength*n*(epsilon_final - epsilon_inflection)/(epsilon_final - epsilon_initiation)
        
    stiffness_inflection = (sigma_inflection - 0.) / (epsilon_inflection - 0.)
    effective_damage = 1 - stiffness_inflection / E1
    stiff_inflection = stiffness6x6(E1, E2, G12, nu12, nu23, d1=effective_damage)
    stretch_inflection = stretch_from_gl_strain(epsilon_inflection)
    return gl_strain_cauchy_stress_from_stretch(stretch_inflection, stiff_inflection)[1][0,0,0]

def displacement_at_inflection(**kwargs):
    args = generic_variable_names_from_dict(**kwargs)
    strength = args['strength']
    E1 = args['E1']
    GX = args['GX']
    n = args['n']
    m = args['m']
    Lc = args['Lc']
    if None in [strength, E1, GX, n, m, Lc]:
        raise KeyError("displacement_at_inflection requires strength, E1, GX, n, m, and Lc to be defined in kwargs.")
    epsilon_initiation = strength/E1
    epsilon_final_A = 2.*GX*m / (strength*n*Lc)  # 'A' represents the first superposed softening law, scaled by n and m
    epsilon_final_B = 2.*GX*(1.-m) / (strength*(1.-n)*Lc)  # 'B' represents the second superposed softening law
    if strength > 0:  # Tension
        epsilon_inflection = max(epsilon_initiation, min(epsilon_final_A, epsilon_final_B))
    else:  # Compression
        epsilon_inflection = min(epsilon_initiation, max(epsilon_final_A, epsilon_final_B))
    stretch_inflection = stretch_from_gl_strain(epsilon_inflection)
    return (stretch_inflection-1)*Lc

def displacement_at_final(**kwargs):
    args = generic_variable_names_from_dict(**kwargs)
    strength = args['strength']
    GX = args['GX']
    n = args['n']
    m = args['m']
    Lc = args['Lc']
    if None in [strength, GX, n, m, Lc]:
        raise KeyError("displacement_at_final requires strength, GX, n, m, and Lc to be defined in kwargs.")
    epsilon_final_A = 2.*GX*m / (strength*n*Lc)  # 'A' represents the first superposed softening law, scaled by n and m
    epsilon_final_B = 2.*GX*(1.-m) / (strength*(1.-n)*Lc)  # 'B' represents the second superposed softening law
    if strength > 0:  # Tension
        epsilon_final = max(epsilon_final_A, epsilon_final_B)
    else:  # Compression
        epsilon_final = min(epsilon_final_A, epsilon_final_B)
    return (stretch_from_gl_strain(epsilon_final)-1)*Lc

def plot_trilinear_law(strength, E1, E2, G12, nu12, nu23, GX, n, m, Lc):
    """
    Note, this only partially accounts for finite strain nonlinearity, and will
    not agree perfectly with Abaqus results when strains are large.
    """
    # Plot input stress-strain curve
    epsilon_initiation = strength/E1
    epsilon_final_A = 2.*GX*m / (strength*n*Lc)  # 'A' represents the first superposed softening law, scaled by n and m
    epsilon_final_B = 2.*GX*(1.-m) / (strength*(1.-n)*Lc)  # 'B' represents the second superposed softening law\
    strength_A = strength*n
    strength_B = strength*(1. - n)
    if strength > 0:  # Tension
        epsilon_inflection = max(epsilon_initiation, min(epsilon_final_A, epsilon_final_B))
        epsilon_final      = max(epsilon_final_A, epsilon_final_B)
    else:
        epsilon_inflection = min(epsilon_initiation, max(epsilon_final_A, epsilon_final_B))
        epsilon_final      = min(epsilon_final_A, epsilon_final_B)
    sigma_inflection = strength*(1. - n)*(epsilon_final - epsilon_inflection)/(epsilon_final - epsilon_initiation)
    plt.plot([0., epsilon_initiation, epsilon_inflection, epsilon_final], [0., strength, sigma_inflection, 0.], 'k', label='Input total')
    plt.plot([0., epsilon_initiation, epsilon_final_A], [0., strength_A, 0.], 'r', label='Input A')
    plt.plot([0., epsilon_initiation, epsilon_final_B], [0., strength_B, 0.], 'b', label='Input B')
    # # Plot stiffness through inflection point
    # stiffness_inflection = (sigma_inflection - 0.) / (epsilon_inflection - 0.)
    # plt.plot([0., epsilon_final], [0., stiffness_inflection*epsilon_final], 'k--', label='Inflection point')
    # # Finite strain nonlinearity curve
    # effective_damage = 1 - stiffness_inflection / E1
    # stiff_inflection = stiffness6x6(E1, E2, G12, nu12, nu23, d1=effective_damage)
    # final_displacement = (stretch_from_gl_strain(epsilon_final)-1)*Lc
    # applied_stretch = np.linspace(1+final_displacement, 1.0, 100)
    # eps11, cauchy = gl_strain_cauchy_stress_from_stretch(applied_stretch, stiff_inflection)
    # plt.plot(eps11, cauchy[0,0,:], 'x-', label='CDM')

    # Plot location of inflection point accounting for finite strain nonlinearity
    params = dict(strength=strength, E1=E1, E2=E2, G12=G12, nu12=nu12, nu23=nu23, GX=GX, n=n, m=m, Lc=Lc)
    cauchy = cauchy_stress_at_inflection(**params)
    plt.scatter(epsilon_inflection, cauchy, label='Inflection point (finite strain)', color='k', zorder=5)


    

    if strength > 0:
        plt.xlim(left=0.)
        plt.ylim(bottom=0.)
    else:
        plt.xlim(right=0.)
        plt.ylim(top=0.)
    plt.xlabel('Strain')
    plt.ylabel('Cauchy stress (MPa)')
    plt.legend()
    plt.show()

def gl_2pk_stress_strain(stretch):
    stiff = stiffness6x6(E1=171420*0.1, E2=9080, G12=5290, nu12=0.32, nu23=0.52, d1=0)
    eps11 = np.zeros_like(stretch)
    cauchy = np.zeros((3,3,len(stretch)))
    for i, s in enumerate(stretch):
        F = np.array([[s, 0., 0.], [0., 1., 0.], [0., 0., 1.]])
        eps = 0.5*(np.matmul(F.T, F) - np.eye(3))
        eps_vec = np.array([eps[0,0], eps[1,1], eps[2,2], 2*eps[0,1], 2*eps[1,2], 2*eps[0,2]])
        pk2_vec = np.matmul(stiff, eps_vec)
        pk2 = np.array([[pk2_vec[0], pk2_vec[3], pk2_vec[5]], [pk2_vec[3], pk2_vec[1], pk2_vec[4]], [pk2_vec[5], pk2_vec[4], pk2_vec[2]]])
        cauchy[:,:,i] = np.matmul(np.matmul(F.T, pk2), F) / np.linalg.det(F)
        eps11[i] = eps[0,0]
    plt.plot(eps11, cauchy[0,0,:], 'k')
    plt.show()


if __name__ == "__main__":

    # Tension example parameters
    # params = dict(
    #     strength = 2326.2,
    #     E1 = 171420.,
    #     E2 = 9080,
    #     G12 = 5290,
    #     nu12 = 0.32,
    #     nu23 = 0.52,
    #     GX = 205.0,
    #     n = 0.63,
    #     m = 0.469,
    #     Lc = 0.5      # Fiber-direction element size
    # )
    # plot_trilinear_law(**params)

    # Compression example parameters
    params = dict(
        strength = -1200.1,
        E1 = 171420.,
        E2 = 9080,
        G12 = 5290,
        nu12 = 0.32,
        nu23 = 0.52,
        GX = 61.0,
        n = 0.614,
        m = 0.432,
        Lc = 0.5      # Fiber-direction element size
    )
    plot_trilinear_law(**params)
