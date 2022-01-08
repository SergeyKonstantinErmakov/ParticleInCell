import numpy as np
import cupy as cp 
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec
from matplotlib.ticker import LinearLocator

import Constants

c_0 = Constants.c_0 # speed of light

N = Constants.N

Nx = Constants.Nx
Ny = Constants.Ny
Nz = Constants.Nz
dx = Constants.dx
dy = Constants.dy
dz = Constants.dz

dt = Constants.dt
t_f = Constants.t_f
print("ok")

for t in range(19000000, 5000000000, 250):
    #-------------------------------------------------------------------------------------
    print("hallo")
    R_pDir = 'C:/Users/serge/OneDrive/Dokumente/Besondere Lernleistung/SimulationCases/QuadrupoleReconnection - RealisticValues/Simulation3_Arrays/R_p_t=' + str(t) + 'dt.npy'
    RPlot_p = cp.asnumpy(cp.load(R_pDir))
    
    V_pDir = 'C:/Users/serge/OneDrive/Dokumente/Besondere Lernleistung/SimulationCases/QuadrupoleReconnection - RealisticValues/Simulation3_Arrays/V_p_t=' + str(t) + 'dt.npy'
    VPlot_p = cp.asnumpy(cp.load(V_pDir))
    _V_ = (VPlot_p[0]**2 + VPlot_p[1]**2 + VPlot_p[2]**2)**0.5/1000

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(RPlot_p[0], RPlot_p[1], RPlot_p[2], facecolors=cm.jet(_V_), s=0.1, marker="v")    
    m = cm.ScalarMappable(cmap=cm.jet)
    m.set_array(_V_)
    cbar = plt.colorbar(m, orientation="vertical")    
    #ax.set_xlim(0, Nx*dx)
    #ax.set_ylim(0, Ny*dy)
    #ax.set_zlim(0, Nz*dz)    
    fig.set_size_inches(13, 7)
    plt.title('Protons at t = ' + str(int(t/dt)) + 'dt', fontsize=15)    
    ax.set_xlabel('$x - Axes$', fontsize=12, rotation=0)
    ax.set_ylabel('$y - Axes$', fontsize=12, rotation=0)
    ax.set_zlabel('$z - Axes$', fontsize=12, rotation=0)
    cbar.ax.set_ylabel('$Durchschnittliche Geschwindigkeit            [km/s]$', rotation=270, fontsize=12, labelpad=17)
    # plt.savefig('C:/Users/serge/OneDrive/Dokumente/Besondere Lernleistung/MagneticReconnectionGPU/simulation1/Protons at t = ' + str(int(t/dt)) +  'dt.png')
    plt.show()
    #------------------------------------------------------------------------------------------------------------
    R_eDir = 'C:/Users/serge/OneDrive/Dokumente/Besondere Lernleistung/SimulationCases/QuadrupoleReconnection - RealisticValues/Simulation3_Arrays/R_e_t=' + str(t) + 'dt.npy'
    RPlot_e = cp.asnumpy(cp.load(R_eDir))
    
    V_eDir = 'C:/Users/serge/OneDrive/Dokumente/Besondere Lernleistung/SimulationCases/QuadrupoleReconnection - RealisticValues/Simulation3_Arrays/V_e_t=' + str(t) + 'dt.npy'
    VPlot_e = cp.asnumpy(cp.load(V_eDir))
    
    _V_ = (VPlot_e[0]**2 + VPlot_e[1]**2 + VPlot_e[2]**2)**0.5/1000
    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(RPlot_e[0], RPlot_e[1], RPlot_e[2], facecolors=cm.jet(_V_), s=0.1)    
    m = cm.ScalarMappable(cmap=cm.jet)
    m.set_array(_V_)
    cbar = plt.colorbar(m, orientation="vertical")    
    #ax.set_xlim(0, Nx*dx)
    #ax.set_ylim(0, Ny*dy)
    #ax.set_zlim(0, Nz*dz)    
    fig.set_size_inches(13, 7)
    plt.title('Electrons at t = ' + str(int(t/dt)) + 'dt', fontsize=15)    
    ax.set_xlabel('$x - Axes$', fontsize=12, rotation=0)
    ax.set_ylabel('$y - Axes$', fontsize=12, rotation=0)
    ax.set_zlabel('$z - Axes$', fontsize=12, rotation=0)
    # plt.savefig('C:/Users/serge/OneDrive/Dokumente/Besondere Lernleistung/MagneticReconnectionGPU/simulation1/Electrons at t = ' + str(int(t/dt)) +  'dt.png')
    cbar.ax.set_ylabel('$Durchschnittliche Geschwindigkeit            [km/s]$', rotation=270, fontsize=12, labelpad=17)
    #-------------------------------------------------------------------------------------------
    BDir = 'C:/Users/serge/OneDrive/Dokumente/Besondere Lernleistung/SimulationCases/QuadrupoleReconnection - RealisticValues/Simulation3_Arrays/B_t=' + str(t) + 'dt.npy'
    B = cp.array(cp.load(BDir))
    
    Y_C = 50
    XX = cp.arange(0, Nx*dx, dx)
    ZZ = cp.arange(0, Nz*dz, dz)
    X_PLOT, Z_PLOT = cp.meshgrid(XX, ZZ)
    B_PLOT = cp.asnumpy((cp.swapaxes(B[0], 0, 1)[Y_C]**2 + cp.swapaxes(B[1], 0, 1)[Y_C]**2 + cp.swapaxes(B[2], 0, 1)[Y_C]**2)**0.5)
    fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
    # Plot the surface.
    surf = ax.plot_surface(cp.asnumpy(X_PLOT), cp.asnumpy(Z_PLOT), B_PLOT, cmap=cm.coolwarm, linewidth=0, antialiased=False)
    # Customize the z axis.
    ax.zaxis.set_major_locator(LinearLocator(10))
    fig.set_size_inches(13, 7)
    # Add a color bar which maps values to colors.
    fig.colorbar(surf, shrink=0.5, aspect=5)
    plt.title('Magnetic Field', fontsize=15)    
    ax.set_xlabel('$x - Axes$', fontsize=12, rotation=0)
    ax.set_ylabel('$z - Axes$', fontsize=12)
    # plt.savefig('C:/Users/serge/OneDrive/Dokumente/Besondere Lernleistung/MagneticReconnectionGPU/simulation1/Magnetic field at t = ' + str(int(t/dt)) +  'dt.png')
    # ax.view_init(60, 35)
    plt.show()    
    #-------------------------------------------------------------------------------------------
    EDir = 'C:/Users/serge/OneDrive/Dokumente/Besondere Lernleistung/SimulationCases/QuadrupoleReconnection - RealisticValues/Simulation3_Arrays/E_t=' + str(t) + 'dt.npy'
    E = cp.array(cp.load(EDir))
    
    Y_C = 50
    X_PLOT, Z_PLOT = cp.meshgrid(XX, ZZ)
    E_PLOT = cp.asnumpy((cp.swapaxes(E[0], 0, 1)[Y_C]**2 + cp.swapaxes(E[1], 0, 1)[Y_C]**2 + cp.swapaxes(E[2], 0, 1)[Y_C]**2)**0.5)
    fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
    # Plot the surface.
    surf = ax.plot_surface(cp.asnumpy(X_PLOT), cp.asnumpy(Z_PLOT), E_PLOT, cmap=cm.coolwarm, linewidth=0, antialiased=False)
    # Customize the z axis.
    ax.zaxis.set_major_locator(LinearLocator(10))
    fig.set_size_inches(13, 7)
    # Add a color bar which maps values to colors.
    fig.colorbar(surf, shrink=0.5, aspect=5)
    plt.title('Electric Field', fontsize=15)    
    ax.set_xlabel('$x - Axes$', fontsize=12, rotation=0)
    ax.set_ylabel('$z - Axes$', fontsize=12)
    # plt.savefig('C:/Users/serge/OneDrive/Dokumente/Besondere Lernleistung/MagneticReconnectionGPU/simulation1/Electric field at t = ' + str(int(t/dt)) +  'dt.png')
    # ax.view_init(60, 35)
    plt.show()    
    #-------------------------------------------------------------------------------------------
    JDir = 'C:/Users/serge/OneDrive/Dokumente/Besondere Lernleistung/SimulationCases/QuadrupoleReconnection - RealisticValues/Simulation3_Arrays/J_t=' + str(t) + 'dt.npy'
    J = cp.array(cp.load(JDir))
    
    Y_C = 50
    X_PLOT, Z_PLOT = cp.meshgrid(XX, ZZ)
    J_PLOT = cp.asnumpy((cp.swapaxes(J[0], 0, 1)[Y_C]**2 + cp.swapaxes(J[1], 0, 1)[Y_C]**2 + cp.swapaxes(J[2], 0, 1)[Y_C]**2)**0.5)
    fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
    # Plot the surface.
    surf = ax.plot_surface(cp.asnumpy(X_PLOT), cp.asnumpy(Z_PLOT), J_PLOT, cmap=cm.coolwarm, linewidth=0, antialiased=False)
    # Customize the z axis.
    ax.zaxis.set_major_locator(LinearLocator(10))
    fig.set_size_inches(13, 7)
    # Add a color bar which maps values to colors.
    fig.colorbar(surf, shrink=0.5, aspect=5)
    plt.title('J Field', fontsize=15)    
    ax.set_xlabel('$x - Axis$', fontsize=12, rotation=0)
    ax.set_ylabel('$z - Axis$', fontsize=12)
    # plt.savefig('C:/Users/serge/OneDrive/Dokumente/Besondere Lernleistung/MagneticReconnectionGPU/simulation1/J field at t = ' + str(int(t/dt)) +  'dt.png')
    # ax.view_init(60, 35)
    plt.show()    
    #-------------------------------------------------------------------------------------------
    DATAX = cp.arange(0, Nx*dx, dx)
    DATAZ = cp.arange(0, Nz*dz, dz)
    DATAX, DATAZ = cp.meshgrid(DATAX, DATAZ)
    B_ex = cp.array(cp.load('Quad.npy'))    
    
    BDir = 'C:/Users/serge/OneDrive/Dokumente/Besondere Lernleistung/SimulationCases/QuadrupoleReconnection - RealisticValues/Simulation3_Arrays/B_t=' + str(t) + 'dt.npy'
    B = cp.array(cp.load(BDir))
    
    BU = cp.swapaxes(cp.swapaxes(B[0], 0, 2)[50], 0, 1) + cp.swapaxes(cp.swapaxes(B_ex[0], 0, 2)[50], 0, 1)
    BV = cp.swapaxes(cp.swapaxes(B[1], 0, 2)[50], 0, 1) + cp.swapaxes(cp.swapaxes(B_ex[1], 0, 2)[50], 0, 1)
    
    B_PLOT = DATAX*3.0
    
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    
    x, y, z = np.meshgrid(np.arange(0, Nx*dx, dx), np.arange(0, Ny*dy, dy), np.arange(0, Nz*dz, dz))
    
    ax.quiver(x, y, z, cp.asnumpy(B[0]), cp.asnumpy(B[1]), cp.asnumpy(B[2]), length=0.1)
    
    plt.show()
    
    
    XX_ = np.array(cp.asnumpy(DATAX))
    ZZ_ = np.array(cp.asnumpy(DATAZ))
    BU_ = np.array(cp.asnumpy(BU))
    BV_ = np.array(cp.asnumpy(BV))
    B_PLOT_ = np.array(cp.asnumpy(B_PLOT))
    # , color=B_PLOT_, linewidth=1, cmap='cool'
    fig = plt.figure()
    fig.set_size_inches(50, 45)
    gs = gridspec.GridSpec(nrows=3, ncols=2, height_ratios=[1, 1, 1])        
    ax1 = fig.add_subplot(gs[0, 1])
    strm = ax1.streamplot(XX_, ZZ_, BU_, BV_, density=[3, 3])
    fig.colorbar(strm.lines)
    ax1.set_title('Magnetic Field Lines')
    #-------------------------------------------------------------------------------------------
    DATAX = cp.arange(0, Nx*dx, dx)
    DATAZ = cp.arange(0, Nz*dz, dz)
    DATAX, DATAZ = cp.meshgrid(DATAX, DATAZ)
    B_ex = cp.array(cp.load('Quad.npy'))    
    BU = cp.swapaxes(B[0], 0, 1)[50] + cp.swapaxes(B_ex[0], 0, 1)[50]
    BV = cp.swapaxes(B[2], 0, 1)[50] + cp.swapaxes(B_ex[2], 0, 1)[50]
    B_PLOT = DATAX*3.0
    
    XX_ = np.array(cp.asnumpy(DATAX))
    ZZ_ = np.array(cp.asnumpy(DATAZ))
    BU_ = np.array(cp.asnumpy(BU))
    BV_ = np.array(cp.asnumpy(BV))
    B_PLOT_ = np.array(cp.asnumpy(B_PLOT))
    # , color=B_PLOT_, linewidth=1, cmap='cool'
    fig = plt.figure()
    fig.set_size_inches(50, 45)
    gs = gridspec.GridSpec(nrows=3, ncols=2, height_ratios=[1, 1, 1])        
    ax1 = fig.add_subplot(gs[0, 1])
    strm = ax1.streamplot(XX_, ZZ_, BU_, BV_, density=[3, 3])
    fig.colorbar(strm.lines)
    ax1.set_title('Magnetic Field Lines')
    #-------------------------------------------------------------------------------------------
