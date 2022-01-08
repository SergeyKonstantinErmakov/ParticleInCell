import cupy as cp

import Fields
import Constants
"""-------------------constants----------------------------------------------------------------------"""
c_0 = Constants.c_0 # speed of light

N = Constants.N

Nx = Constants.Nx
Ny = Constants.Ny
Nz = Constants.Nz
dx = Constants.dx
dy = Constants.dy
dz = Constants.dz

dt = Constants.dt
"""-------------------constants----------------------------------------------------------------------"""

_1_ = cp.ones(N)



# for relativistic velocities
def Y_Factor(beta):
    return cp.reciprocal((cp.array(_1_ - beta**2))**0.5)




def Particle(q, m, n_0, R, V, E_0, B_0):

    E = E_0 + Fields.E_external(R)
    B = B_0 + Fields.B_external(R)
   
    #move through phase space
    _V_ = (V[0]**2 + V[1]**2 + V[2]**2)**0.5
    
    U_x = Y_Factor(_V_/c_0)*V[0]
    U_y = Y_Factor(_V_/c_0)*V[1]
    U_z = Y_Factor(_V_/c_0)*V[2]
    
    U_x_s = U_x + q*dt/m*(E[0] + 0.5*V[1]*B[2] - 0.5*V[2]*B[1]) 
    U_y_s = U_y + q*dt/m*(E[1] + 0.5*V[2]*B[0] - 0.5*V[0]*B[2]) 
    U_z_s = U_z + q*dt/m*(E[2] + 0.5*V[0]*B[1] - 0.5*V[1]*B[0])
 
    tau = q*dt/(2*m)*B
    w = (U_x_s*tau[0] + U_y_s*tau[1] + U_z_s*tau[2])/c_0
    gamma_ = (_1_ + (U_x_s**2 + U_y_s**2 + U_z_s**2)**2/c_0**2)**0.5
    sig = gamma_**2 - tau[0]**2 - tau[1]**2 - tau[2]**2
    gamma = (sig + (sig**2 + 4*(tau[0]**2 + tau[1]**2 + tau[2]**2 + w**2))**0.5)**0.5/2
    t_vector = tau/gamma
    
    U_x = (U_x_s + (U_x_s*t_vector[0] + U_y_s*t_vector[1] + U_z_s*t_vector[2])*t_vector[0] + U_y_s*t_vector[2] - U_z_s*t_vector[1])/(1 + t_vector[0]**2 + t_vector[1]**2 + t_vector[2]**2)
    U_y = (U_y_s + (U_x_s*t_vector[0] + U_y_s*t_vector[1] + U_z_s*t_vector[2])*t_vector[1] + U_z_s*t_vector[0] - U_x_s*t_vector[2])/(1 + t_vector[0]**2 + t_vector[1]**2 + t_vector[2]**2)
    U_z = (U_z_s + (U_x_s*t_vector[0] + U_y_s*t_vector[1] + U_z_s*t_vector[2])*t_vector[2] + U_x_s*t_vector[1] - U_y_s*t_vector[0])/(1 + t_vector[0]**2 + t_vector[1]**2 + t_vector[2]**2)

    V[0] = cp.reciprocal(Y_Factor(_V_/c_0))*U_x 
    V[1] = cp.reciprocal(Y_Factor(_V_/c_0))*U_y
    V[2] = cp.reciprocal(Y_Factor(_V_/c_0))*U_z
    
    R[0] += V[0]*dt/2 
    R[1] += V[1]*dt/2
    R[2] += V[2]*dt/2
    
    
    return [R, V]
    










