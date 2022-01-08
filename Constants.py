c_0 = 2.99*10**8 # speed of light
kB = 1.38*10**-23 # boltzmann constant
e_0 = 8.854*10**-12 # electric field constant
u_0 = 1.256*10**-6 # magnetic field constant


n_e = 10**4 # number of negative particles per macroparticle
n_p = 10**4 # number of positive particles per macroparticle
e_ = 1.602*10**-19 # elemtary charge
m_e = 9.109*10**-31 # mass of negative particles 
m_p = 6.646*10**-27 # mass of positive particles 


Nx = 100 # number of cells in x
Ny = 100 # number of cells in x
Nz = 100 # number of cells in x
dx = 2.5*10**-3 # stepsize x
dy = 2.5*10**-3 # stepsize y
dz = 2.5*10**-3 # stepsize z


dt = 10**-16 # stepsize t
t_i = 0 # initial t
t_f = 1000 # final t


T_i = 293 # ion temperature
T_e = 293 # electron temperature


Nx_ = 50 # reconnection box in x direction
Ny_ = 50 # reconnection box in y direction 
Nz_ = 50 # reconnection box in z direction
L = 0.25*Nz_*dz 
L_x = Nx_*dx
L_z = Nz_*dz
B_zero = 1


N = 1000000 # Number of particles
dx_p = dx # Particle destribution size in x
dy_p = dy # Particle destribution size in x
dz_p = dz # Particle destribution size in x