import cupy as cp

import Constants
import Fields
"""-------------------constants----------------------------------------------------------------------"""
pi = cp.pi

N = Constants.N # Number of particles
dx_p = Constants.dx_p # Particle destribution size in x
dy_p = Constants.dy_p # Particle destribution size in x
dz_p = Constants.dz_p # Particle destribution size in x

Nx = Constants.Nx # number of cells in x
Ny = Constants.Ny # number of cells in x
Nz = Constants.Nz # number of cells in x
dx = Constants.dx # stepsize x
dy = Constants.dy # stepsize y
dz = Constants.dz # stepsize z

dt = Constants.dt

n_e = Constants.n_e # number of negative particles per macroparticle
n_p = Constants.n_p # number of positive particles per macroparticle

NP = Fields.NP
NE = Fields.NE

e_ = Constants.e_ # elemtary charge

e_0 = Constants.e_0 # electric field constant
u_0 = Constants.u_0 # magnetic field constant
"""-------------------constants----------------------------------------------------------------------"""

_1_ = cp.ones(N)





# create fft frequencies for FDTD solver
K1 = 2*pi*cp.transpose(cp.array([[cp.fft.fftfreq(Nx, d = dx)]])) # note that the methods of transposition are applied in order to get correct dimensions in the FTDT solver
K2 = 2*pi*cp.transpose(cp.array([cp.fft.fftfreq(Ny, d = dy)])) # (the dimensions must be adjusted as we are performing matrix multiplication)
K3 = 2*pi*cp.array(cp.fft.fftfreq(Nz, d = dz))
K = [K1, K2, K3]







def MaxwellSolver(R_e, V_e, R_p, V_p, E, B, J):
    RINT_e = cp.array([cp.floor((cp.array(R_e[0])/dx).astype(int)).astype(int), cp.floor((cp.array(R_e[1])/dy).astype(int)).astype(int), cp.floor((cp.array(R_e[2])/dz).astype(int)).astype(int)])
    RINT_0_e = cp.array([cp.floor(((cp.array(R_e[0]) - dx_p)/dx).astype(int)).astype(int), cp.floor(((cp.array(R_e[1]) - dy_p)/dy).astype(int)).astype(int), cp.floor(((cp.array(R_e[2]) - dz_p)/dz).astype(int)).astype(int)])
    RINT_1_e = cp.array([cp.floor(((cp.array(R_e[0]) + dx_p)/dx).astype(int)).astype(int), cp.floor(((cp.array(R_e[1]) + dy_p)/dy).astype(int)).astype(int), cp.floor(((cp.array(R_e[2]) + dz_p)/dz).astype(int)).astype(int)])
    RINT_E = [RINT_0_e, RINT_1_e]
    
    RINT_p = cp.array([cp.floor((cp.array(R_p[0])/dx).astype(int)).astype(int), cp.floor((cp.array(R_p[1])/dy).astype(int)).astype(int), cp.floor((cp.array(R_p[2])/dz).astype(int)).astype(int)])
    RINT_0_p = cp.array([cp.floor(((cp.array(R_p[0]) - dx_p)/dx).astype(int)).astype(int), cp.floor(((cp.array(R_p[1]) - dy_p)/dy).astype(int)).astype(int), cp.floor(((cp.array(R_p[2]) - dz_p)/dz).astype(int)).astype(int)])
    RINT_1_p = cp.array([cp.floor(((cp.array(R_p[0]) + dx_p)/dx).astype(int)).astype(int), cp.floor(((cp.array(R_p[1]) + dy_p)/dy).astype(int)).astype(int), cp.floor(((cp.array(R_p[2]) + dz_p)/dz).astype(int)).astype(int)])
    RINT_P = [RINT_0_p, RINT_1_p]
    
    # divide each element (0,1,2) by dx, dy or dz...
    OV0_e = 3*(cp.array([_1_*dx_p/dx, _1_*dy_p/dy, _1_*dz_p/dz]) - cp.array([cp.array(R_e[0])/dx + _1_, cp.array(R_e[1])/dy + _1_, cp.array(R_e[2])/dz + _1_]) - RINT_1_e)
    OV1_e = OV0_e
    
    OV0_p = 3*(cp.array([_1_*dx_p/dx, _1_*dy_p/dy, _1_*dz_p/dz]) - cp.array([cp.array(R_p[0])/dx + _1_, cp.array(R_p[1])/dy + _1_, cp.array(R_p[2])/dz + _1_]) - RINT_1_p)
    OV1_p = OV0_p
    
    #-------------------locating particles-----------------------------------------------------------
    Weights_e = -e_*n_e*NE*cp.array(V_e)/(dx*dy*dz)
    Weights_p = e_*n_p*NP*cp.array(V_p)/(dx*dy*dz)
    #-------------------for J calculation------------------------------------------------------------
    
    
    # for all three dimensions
    SPLNe = []
    SPLNp = []
    
    SPLN_e_0 = [cp.ones(N), cp.ones(N), cp.ones(N)]
    SPLN_e_1 = [cp.ones(N), cp.ones(N), cp.ones(N)]
    SPLN_p_0 = [cp.ones(N), cp.ones(N), cp.ones(N)]
    SPLN_p_1 = [cp.ones(N), cp.ones(N), cp.ones(N)]
    
    
    for d in range(3):
        
        # back----------
        # if 0 <= t and t < 1:
        T1 = cp.where(0 <= OV0_e[d])[0]
        T2 = cp.where(OV0_e[d] < 1)[0]
        T = T1[cp.in1d(T1, T2)]
        SPLN_e_0[d][T] *= (1/6*OV0_e[d][T]**3).astype(float)        
        # elif 1 <= t and t < 2:
        T1 = cp.where(1 <= OV0_e[d])[0]
        T2 = cp.where(OV0_e[d] < 2)[0]
        T = T1[cp.in1d(T1, T2)]
        SPLN_e_0[d][T] *= (0.5*_1_[T] - 1/3*OV0_e[d][T]**3 + 3/2*OV0_e[d][T]**2 - 3/2*OV0_e[d][T]).astype(float)
        # elif 2 <= t:
        T = cp.where(2 <= OV0_e[d])[0]
        SPLN_e_0[d][T] *= (1/6*OV0_e[d][T]**3 - 3/2*OV0_e[d][T]**2 + 9/2*OV0_e[d][T] - 7/2*_1_[T]).astype(float)        
        
        # front-----------
        # if 2 < t:
        T = cp.where(2 < OV1_e[d])
        SPLN_e_1[d][T] *= (4.5*_1_[T] - 1/6*OV1_e[d][T]**3 + 3/2*OV1_e[d][T]**2 - 4.5*OV1_e[d][T]).astype(float)
        # elif 1 < t and t <= 2:
        T1 = cp.where(1 < OV1_e[d])[0]
        T2 = cp.where(OV1_e[d] <= 2)[0]
        T = T1[cp.in1d(T1, T2)]
        SPLN_e_1[d][T] *= (0.5*_1_[T] + 1/3*OV1_e[d][T]**3 - 3/2*OV1_e[d][T]**2 + 3/2*OV1_e[d][T]).astype(float)
        # elif 0 <= t and t <= 1:
        T1 = cp.where(0 <= OV1_e[d])[0]
        T2 = cp.where(OV1_e[d] <= 1)[0]
        T = T1[cp.in1d(T1, T2)]
        SPLN_e_1[d][T] *= (1*_1_[T] - 1/6*OV1_e[d][T]**3).astype(float)    
        
        # back----------
        # if 0 <= t and t < 1:
        T1 = cp.where(0 <= OV0_p[d])[0]
        T2 = cp.where(OV0_p[d] < 1)[0]
        T = T1[cp.in1d(T1, T2)]
        SPLN_p_0[d][T] *= (1/6*OV0_p[d][T]**3).astype(float)
        # elif 1 <= t and t < 2:
        T1 = cp.where(1 <= OV0_p[d])[0]
        T2 = cp.where(OV0_p[d] < 2)[0]
        T = T1[cp.in1d(T1, T2)]
        SPLN_p_0[d][T] = (0.5*_1_[T] - 1/3*OV0_p[d][T]**3 + 3/2*OV0_p[d][T]**2 - 3/2*OV0_p[d][T]).astype(float)
        # elif 2 <= t:
        T = cp.where(2 < OV0_p[d])[0]
        SPLN_p_0[d][T] *= (1/6*OV0_p[d][T]**3 - 3/2*OV0_p[d][T]**2 + 9/2*OV0_p[d][T] - 7/2*_1_[T]).astype(float)
        
        # front-----------
        # if 2 < t:
        T = cp.where(2 < OV1_p[d])[0]
        SPLN_p_1[d][T] *= (4.5*_1_[T] - 1/6*OV1_p[d][T]**3 + 3/2*OV1_p[d][T]**2 - 4.5*OV1_p[d][T]).astype(float)
        # elif 1 < t and t <= 2:
        T1 = cp.where(1 < OV1_p[d])[0]
        T2 = cp.where(OV1_p[d] <= 2)[0]
        T = T1[cp.in1d(T1, T2)]
        SPLN_p_1[d][T] *= (0.5*_1_[T] + 1/3*OV1_p[d][T]**3 - 3/2*OV1_p[d][T]**2 + 3/2*OV1_p[d][T]).astype(float)
        # elif 0 <= t and t <= 1:
        T1 = cp.where(0 <= OV1_p[d])[0]
        T2 = cp.where(OV1_p[d] <= 1)[0]
        T = T1[cp.in1d(T1, T2)]
        SPLN_p_1[d][T] = (1*_1_[T] - 1/6*OV1_e[d][T]**3).astype(float)    
        
        SPLNe.append([SPLN_e_0[d], SPLN_e_1[d]])
        SPLNp.append([SPLN_p_0[d], SPLN_p_1[d]])
        
        
    
    # all overlap combinations of back and front
    for kappax in range(2):
        for kappay in range(2):
            for kappaz in range(2):
                # electron contribution
                J[0][RINT_E[kappax][0], RINT_E[kappay][1], RINT_E[kappaz][2]] += cp.array(SPLNe[0][kappax]*SPLNe[1][kappay]*SPLNe[2][kappaz]*Weights_e[0])
                # proton contribution
                J[0][RINT_P[kappax][0], RINT_P[kappay][1], RINT_P[kappaz][2]] += cp.array(SPLNp[0][kappax]*SPLNp[1][kappay]*SPLNp[2][kappaz]*Weights_p[0])
                
                J[1][RINT_E[kappax][0], RINT_E[kappay][1], RINT_E[kappaz][2]] += cp.array(SPLNe[0][kappax]*SPLNe[1][kappay]*SPLNe[2][kappaz]*Weights_e[1])
                J[1][RINT_P[kappax][0], RINT_P[kappay][1], RINT_P[kappaz][2]] += cp.array(SPLNp[0][kappax]*SPLNp[1][kappay]*SPLNp[2][kappaz]*Weights_p[1])
                
                J[2][RINT_E[kappax][0], RINT_E[kappay][1], RINT_E[kappaz][2]] += cp.array(SPLNe[0][kappax]*SPLNe[1][kappay]*SPLNe[2][kappaz]*Weights_e[2])
                J[2][RINT_P[kappax][0], RINT_P[kappay][1], RINT_P[kappaz][2]] += cp.array(SPLNp[0][kappax]*SPLNp[1][kappay]*SPLNp[2][kappaz]*Weights_p[2])
                
        
    
    # FDTD method for solving Maxwell's Equations
    FE = [cp.fft.fftn(E[0]), cp.fft.fftn(E[1]), cp.fft.fftn(E[2])]
    FB = [cp.fft.fftn(B[0]), cp.fft.fftn(B[1]), cp.fft.fftn(B[2])]                
    FJ = [cp.fft.fftn(J[0]), cp.fft.fftn(J[1]), cp.fft.fftn(J[2])]
    
    FE[0] += 1/e_0*(1j*dt*(K[1]*FB[2] - K[2]*FB[1]) - dt*FJ[0])    
    FE[1] += 1/e_0*(1j*dt*(K[2]*FB[0] - K[0]*FB[2]) - dt*FJ[1])    
    FE[2] += 1/e_0*(1j*dt*(K[0]*FB[1] - K[1]*FB[0]) - dt*FJ[2])    
    
    FB[0] -= 1/u_0*(1j*dt*(K[1]*FE[2] - K[2]*FE[1]))
    FB[1] -= 1/u_0*(1j*dt*(K[2]*FE[0] - K[0]*FE[2]))
    FB[2] -= 1/u_0*(1j*dt*(K[0]*FE[1] - K[1]*FE[0]))
                
    E = [cp.fft.ifftn(FE[0]).real, cp.fft.ifftn(FE[1]).real, cp.fft.ifftn(FE[2]).real]
    B = [cp.fft.ifftn(FB[0]).real, cp.fft.ifftn(FB[1]).real, cp.fft.ifftn(FB[2]).real]                
    J = [cp.fft.ifftn(FJ[0]).real, cp.fft.ifftn(FJ[1]).real, cp.fft.ifftn(FJ[2]).real]
    
    
    # allign as suitable array for particles
    E_e = cp.array([E[0][RINT_e[0], RINT_e[1], RINT_e[2]], E[1][RINT_e[0], RINT_e[1], RINT_e[2]], E[2][RINT_e[0], RINT_e[1], RINT_e[2]]])
    B_e = cp.array([B[0][RINT_e[0], RINT_e[1], RINT_e[2]], B[1][RINT_e[0], RINT_e[1], RINT_e[2]], B[2][RINT_e[0], RINT_e[1], RINT_e[2]]])
    E_p = cp.array([E[0][RINT_p[0], RINT_p[1], RINT_p[2]], E[1][RINT_p[0], RINT_p[1], RINT_p[2]], E[2][RINT_p[0], RINT_p[1], RINT_p[2]]])
    B_p = cp.array([B[0][RINT_p[0], RINT_p[1], RINT_p[2]], B[1][RINT_p[0], RINT_p[1], RINT_p[2]], B[2][RINT_p[0], RINT_p[1], RINT_p[2]]])
    
    return [E, B, E_e, B_e, E_p, B_p]

    
