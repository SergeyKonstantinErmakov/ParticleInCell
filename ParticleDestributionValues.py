import cupy as cp

import Constants
from scipy.stats import maxwell

Nx = Constants.Nx # number of cells in x
Ny = Constants.Ny # number of cells in x
Nz = Constants.Nz # number of cells in x
dx = Constants.dx # stepsize x
dy = Constants.dy # stepsize y
dz = Constants.dz # stepsize z

N = Constants.N # Number of particles

T_i = Constants.T_i # ion temperature
T_e = Constants.T_e # electron temperature

kB = Constants.kB

m_e = Constants.m_e # mass of negative particles 
m_p = Constants.m_p # mass of positive particles 

R_e = cp.array([cp.random.normal(0.5*Nx*dx, 1/10*Nx*dx, size=N), cp.random.normal(0.5*Ny*dy, 1/10*Ny*dy, size=N), cp.random.normal(0.5*Nz*dz, 1/10*Nz*dz, size=N)]) # position distribution is a gaussian
V_e = cp.array([3**-0.5*maxwell.rvs((kB*T_e/m_e)**0.5, size=N), 3**-0.5*maxwell.rvs((kB*T_e/m_e)**0.5, size=N), 3**-0.5*maxwell.rvs((kB*T_e/m_e)**0.5, size=N)])*10**-4 # velocity maxwellian
R_p = R_e.copy()
V_p = cp.array([3**-0.5*maxwell.rvs((kB*T_i/m_p)**0.5, size=N), 3**-0.5*maxwell.rvs((kB*T_i/m_p)**0.5, size=N), 3**-0.5*maxwell.rvs((kB*T_i/m_p)**0.5, size=N)])*10**-4


cp.save('R_e', R_e)
cp.save('V_e', V_e)
cp.save('R_p', R_p)
cp.save('V_p', V_p)
