import concurrent.futures
import cupy as cp
import numpy as np

import ParticlePusher
import MaxwellSolver
import Fields
import Constants
"""-------------------constants----------------------------------------------------------------------"""
Nx = Constants.Nx # number of cells in x
Ny = Constants.Ny # number of cells in x
Nz = Constants.Nz # number of cells in x
dx = Constants.dx # stepsize x
dy = Constants.dy # stepsize y
dz = Constants.dz # stepsize z

t_i = Constants.t_i # initial t
t_f = Constants.t_f # final t
dt = Constants.dt # stepsize t

n_e = Constants.n_e # number of negative particles per macroparticle
n_p = Constants.n_p # number of positive particles per macroparticle
e_ = Constants.e_ # elementary charge
m_e = Constants.m_e # mass of negative particles 
m_p = Constants.m_p # mass of positive particles 
"""-------------------constants----------------------------------------------------------------------"""



# create a grid 
XX = cp.arange(0, Nx*dx, dx)
YY = cp.arange(0, Ny*dy, dy)
ZZ = cp.arange(0, Nz*dz, dz)

X, Y, Z = cp.meshgrid(XX, YY, ZZ)



R_e = Fields.R_e
V_e = Fields.V_e
R_p = Fields.R_p
V_p = Fields.V_p

E = Fields.EInitial(X, Y, Z)
B = Fields.BInitial(X, Y, Z)
J = Fields.JInitial(X, Y, Z)





t = t_i
while t <= t_f:
    print("t =",t)
    """plotting section"""
    if int(t/dt)%10**5==0:
        BDir = 'C:/Users/serge/OneDrive/Dokumente/Besondere Lernleistung/SimulationCases/QuadrupoleReconnection - RealisticValues/Simulation3_Arrays/B_t=' + str(int(t/dt)/1000) + 'dt'
        BNP = cp.asnumpy(cp.array(B))
        cp.save(BDir, BNP)
        
        EDir = 'C:/Users/serge/OneDrive/Dokumente/Besondere Lernleistung/SimulationCases/QuadrupoleReconnection - RealisticValues/Simulation3_Arrays/E_t=' + str(int(t/dt)/1000) + 'dt'
        ENP = cp.asnumpy(cp.array(E))
        cp.save(EDir, ENP)
        
        JDir = 'C:/Users/serge/OneDrive/Dokumente/Besondere Lernleistung/SimulationCases/QuadrupoleReconnection - RealisticValues/Simulation3_Arrays/J_t=' + str(int(t/dt)/1000) + 'dt'
        JNP = cp.asnumpy(cp.array(J))
        cp.save(JDir, JNP)
        
        R_eDir = 'C:/Users/serge/OneDrive/Dokumente/Besondere Lernleistung/SimulationCases/QuadrupoleReconnection - RealisticValues/Simulation3_Arrays/R_e_t=' + str(int(t/dt)/1000) + 'dt'
        R_eNP = cp.asnumpy(R_e)
        cp.save(R_eDir, R_eNP)
        
        V_eDir = 'C:/Users/serge/OneDrive/Dokumente/Besondere Lernleistung/SimulationCases/QuadrupoleReconnection - RealisticValues/Simulation3_Arrays/V_e_t=' + str(int(t/dt)/1000) + 'dt'
        V_eNP = cp.asnumpy(V_e)
        cp.save(V_eDir, V_eNP)
        
        R_pDir = 'C:/Users/serge/OneDrive/Dokumente/Besondere Lernleistung/SimulationCases/QuadrupoleReconnection - RealisticValues/Simulation3_Arrays/R_p_t=' + str(int(t/dt)/1000) + 'dt'
        R_pNP = cp.asnumpy(R_p)
        cp.save(R_pDir, R_pNP)
        print(R_p)
        V_pDir = 'C:/Users/serge/OneDrive/Dokumente/Besondere Lernleistung/SimulationCases/QuadrupoleReconnection - RealisticValues/Simulation3_Arrays/V_p_t=' + str(int(t/dt)/1000) + 'dt'
        V_pNP = cp.asnumpy(V_p)
        cp.save(V_pDir, V_pNP)
        
    
    # leapfrog pusher
    push = 0
    while push < 2:
        EM = MaxwellSolver.MaxwellSolver(R_e, V_e, R_p, V_p, E, B, J)
        E = EM[0]
        B = EM[1]
        E_e = EM[2]
        B_e = EM[3]
        E_p = EM[4]
        B_p = EM[5]
        
        # threading
        with concurrent.futures.ThreadPoolExecutor() as executor:                            
            electron = executor.submit(ParticlePusher.Particle, -e_, m_e, n_e, R_e, V_e, E_e, B_e)
            proton = executor.submit(ParticlePusher.Particle, e_, m_p, n_p, R_p, V_p, E_p, B_p)
            
            R_e = electron.result()[0]
            V_e = electron.result()[1]
            R_p = proton.result()[0]
            V_p = proton.result()[1]
            
        push += 1
    
    
    
    t += dt

















