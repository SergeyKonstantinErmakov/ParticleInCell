import cupy as cp

import Constants
"""-------------------constants----------------------------------------------------------------------"""
pi = cp.pi

Nx = Constants.Nx # number of cells in x
Ny = Constants.Ny # number of cells in x
Nz = Constants.Nz # number of cells in x
dx = Constants.dx # stepsize x
dy = Constants.dy # stepsize y
dz = Constants.dz # stepsize z

N = Constants.N # Number of particles

Nx_ = Constants.Nx_ # reconnection box in x direction
Ny_ = Constants.Ny_ # reconnection box in y direction 
Nz_ = Constants.Nz_ # reconnection box in z direction
L = Constants.L
L_x = Constants.L_x
L_z = Constants.L_z
B_zero = Constants.B_zero

T_i = Constants.T_i # ion temperature
T_e = Constants.T_e # electron temperature

u_0 = Constants.u_0
kB = Constants.kB

m_e = Constants.m_e # mass of negative particles 
m_p = Constants.m_p # mass of positive particles 
"""-------------------constants----------------------------------------------------------------------"""

_0_ = cp.zeros(N)
_1_ = cp.ones(N)


# define initial particle-/velocity destribution
R_p = cp.array(cp.load('C:/Users/serge/OneDrive/Dokumente/Besondere Lernleistung/SimulationCases/QuadrupoleReconnection - RealisticValues/Simulation3_Arrays/R_p_t=t + 1900000010^-3dt.npy'))    
V_p = cp.array(cp.load('C:/Users/serge/OneDrive/Dokumente/Besondere Lernleistung/SimulationCases/QuadrupoleReconnection - RealisticValues/Simulation3_Arrays/V_p_t=t + 1900000010^-3dt.npy'))    
R_e = cp.array(cp.load('C:/Users/serge/OneDrive/Dokumente/Besondere Lernleistung/SimulationCases/QuadrupoleReconnection - RealisticValues/Simulation3_Arrays/R_e_t=t + 1900000010^-3dt.npy'))    
V_e = cp.array(cp.load('C:/Users/serge/OneDrive/Dokumente/Besondere Lernleistung/SimulationCases/QuadrupoleReconnection - RealisticValues/Simulation3_Arrays/V_e_t=t + 1900000010^-3dt.npy'))    

# initial destribution function! - particle weighing
NP = cp.array(cp.load())
NE = cp.array(cp.load())


# define initial fields
def EInitial(X, Y, Z):
    return cp.array(cp.load('C:/Users/serge/OneDrive/Dokumente/Besondere Lernleistung/SimulationCases/QuadrupoleReconnection - RealisticValues/Simulation3_Arrays/E_t=t + 1900000010^-3dt.npy'))

B_ex = cp.array(cp.load('Quad.npy'))    
def BInitial(X, Y, Z):
    return cp.array(cp.load('C:/Users/serge/OneDrive/Dokumente/Besondere Lernleistung/SimulationCases/QuadrupoleReconnection - RealisticValues/Simulation3_Arrays/B_t=t + 1900000010^-3dt.npy'))
    
def JInitial(X, Y, Z):
    return cp.array(cp.load('C:/Users/serge/OneDrive/Dokumente/Besondere Lernleistung/SimulationCases/QuadrupoleReconnection - RealisticValues/Simulation3_Arrays/J_t=t + 1900000010^-3dt.npy'))
    

# define external Fields
def E_external(R):
    E = cp.array([_0_, _0_, _0_]) # external field values saved as a vector in order to immediately use with particle array
    return E


def B_external(R):
    IntX = cp.floor((cp.array(R[0])/dx).astype(int)).astype(int)
    IntY = cp.floor((cp.array(R[1])/dy).astype(int)).astype(int)
    IntZ = cp.floor((cp.array(R[2])/dz).astype(int)).astype(int)
    B = cp.array([10000.*B_ex[0][IntX, IntY, IntZ], 10000.*B_ex[1][IntX, IntY, IntZ], 10000.*B_ex[2][IntX, IntY, IntZ]])
    return B



