import cupy as cp
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import Constants 


u_0 = Constants.u_0


Nx = Constants.Nx # number of cells in x
Ny = Constants.Ny # number of cells in x
Nz = Constants.Nz # number of cells in x
dx = Constants.dx # stepsize x
dy = Constants.dy # stepsize y
dz = Constants.dz # stepsize z


dphi = 2*cp.pi*10**-3
PHI = cp.arange(0, 2*cp.pi, dphi)



def B(x, y, z, d, R0, I, alpha):
    """--------------------------------------------------Ring Nr. 1,2----------------------------------------------------"""
    B = cp.array([X*0., Y*0., Z*0.])
    
    # summation of wire elements
    for phi in PHI:        
       
        L_x_notRotated = 0
        L_y_notRotated = R0*(cp.cos(phi + dphi) - cp.cos(phi - dphi))*R0*cp.sqrt(2*(1 - cp.cos(dphi)))
        L_z_notRotated = R0*(cp.sin(phi + dphi) - cp.sin(phi - dphi))*R0*cp.sqrt(2*(1 - cp.cos(dphi))) 
        
        # rotation by alpha around z axes
        L_x = L_x_notRotated*cp.cos(alpha) - L_y_notRotated*cp.sin(alpha)
        L_y = L_x_notRotated*cp.sin(alpha) + L_y_notRotated*cp.cos(alpha)
        L_z = L_z_notRotated
        
        r_x_notRotated = d 
        r_y_notRotated = R0*cp.cos(phi) 
        r_z_notRotated = R0*cp.sin(phi)  
        
        # rotation by alpha around z axes as well as subtraction of particle position
        r_x = r_x_notRotated*cp.cos(alpha) - r_y_notRotated*cp.sin(alpha) - x
        r_y = r_x_notRotated*cp.sin(alpha) + r_y_notRotated*cp.cos(alpha) - y
        r_z = r_z_notRotated - z
        
        B_x_A = u_0*I/(4*cp.pi)*(L_y*r_z - L_z*r_y)/((L_x**2 + L_y**2 + L_z**2)*(r_x**2 + r_y**2 + r_z**2) - (r_x*L_x + r_y*L_y + r_z*L_z)**2)
        B_x_B = (L_x**2 + L_y**2 + L_z**2 - r_x*L_x - r_y*L_y - r_z*L_z)/cp.sqrt(r_x**2 + r_y**2 + r_z**2 + 2*(r_x*L_x + r_y*L_y + r_z*L_z) + L_x**2 + L_y**2 + L_z**2) + (r_x*L_x + r_y*L_y + r_z*L_z)/cp.sqrt(r_x**2 + r_y**2 + r_z**2)
        B_x = B_x_A*B_x_B
        
        B_y_A = u_0*I/(4*cp.pi)*(L_z*r_x - L_x*r_z)/((L_x**2 + L_y**2 + L_z**2)*(r_x**2 + r_y**2 + r_z**2) - (r_x*L_x + r_y*L_y + r_z*L_z)**2)
        B_y_B = (L_x**2 + L_y**2 + L_z**2 - r_x*L_x - r_y*L_y - r_z*L_z)/cp.sqrt(r_x**2 + r_y**2 + r_z**2 + 2*(r_x*L_x + r_y*L_y + r_z*L_z) + L_x**2 + L_y**2 + L_z**2) + (r_x*L_x + r_y*L_y + r_z*L_z)/cp.sqrt(r_x**2 + r_y**2 + r_z**2)
        B_y = B_y_A*B_y_B
        
        B_z_A = u_0*I/(4*cp.pi)*(L_x*r_y - L_y*r_x)/((L_x**2 + L_y**2 + L_z**2)*(r_x**2 + r_y**2 + r_z**2) - (r_x*L_x + r_y*L_y + r_z*L_z)**2)
        B_z_B = (L_x**2 + L_y**2 + L_z**2 - r_x*L_x - r_y*L_y - r_z*L_z)/cp.sqrt(r_x**2 + r_y**2 + r_z**2 + 2*(r_x*L_x + r_y*L_y + r_z*L_z) + L_x**2 + L_y**2 + L_z**2) + (r_x*L_x + r_y*L_y + r_z*L_z)/cp.sqrt(r_x**2 + r_y**2 + r_z**2)
        B_z = B_z_A*B_z_B
        
        B += cp.array([B_x, B_y, B_z])
        
    return B




def Borth(x, y, z, d, R0, I, alpha):
    """--------------------------------------------------Ring Nr. 3,4----------------------------------------------------"""
    B = cp.array([X*0., Y*0., Z*0.])
    
    # summation of wire elements
    for phi in PHI:        
        
        L_x_notRotated = R0*(cp.cos(phi + dphi) - cp.cos(phi - dphi))*R0*cp.sqrt(2*(1 - cp.cos(dphi)))
        L_y_notRotated = 0
        L_z_notRotated = R0*(cp.sin(phi + dphi) - cp.sin(phi - dphi))*R0*cp.sqrt(2*(1 - cp.cos(dphi))) 
        
        # rotation by alpha around z axes
        L_x = L_x_notRotated*cp.cos(alpha) - L_y_notRotated*cp.sin(alpha)
        L_y = L_x_notRotated*cp.sin(alpha) + L_y_notRotated*cp.cos(alpha)
        L_z = L_z_notRotated
        
        r_x_notRotated = R0*cp.cos(phi) 
        r_y_notRotated = d 
        r_z_notRotated = R0*cp.sin(phi) 
        
        # rotation by alpha around z axes as well as subtraction of particle position
        r_x = r_x_notRotated*cp.cos(alpha) - r_y_notRotated*cp.sin(alpha) - x
        r_y = r_x_notRotated*cp.sin(alpha) + r_y_notRotated*cp.cos(alpha) - y
        r_z = r_z_notRotated - z
        
        B_x_A = u_0*I/(4*cp.pi)*(L_y*r_z - L_z*r_y)/((L_x**2 + L_y**2 + L_z**2)*(r_x**2 + r_y**2 + r_z**2) - (r_x*L_x + r_y*L_y + r_z*L_z)**2)
        B_x_B = (L_x**2 + L_y**2 + L_z**2 - r_x*L_x - r_y*L_y - r_z*L_z)/cp.sqrt(r_x**2 + r_y**2 + r_z**2 + 2*(r_x*L_x + r_y*L_y + r_z*L_z) + L_x**2 + L_y**2 + L_z**2) + (r_x*L_x + r_y*L_y + r_z*L_z)/cp.sqrt(r_x**2 + r_y**2 + r_z**2)
        B_x = B_x_A*B_x_B
    
        
        B_y_A = u_0*I/(4*cp.pi)*(L_z*r_x - L_x*r_z)/((L_x**2 + L_y**2 + L_z**2)*(r_x**2 + r_y**2 + r_z**2) - (r_x*L_x + r_y*L_y + r_z*L_z)**2)
        B_y_B = (L_x**2 + L_y**2 + L_z**2 - r_x*L_x - r_y*L_y - r_z*L_z)/cp.sqrt(r_x**2 + r_y**2 + r_z**2 + 2*(r_x*L_x + r_y*L_y + r_z*L_z) + L_x**2 + L_y**2 + L_z**2) + (r_x*L_x + r_y*L_y + r_z*L_z)/cp.sqrt(r_x**2 + r_y**2 + r_z**2)
        B_y = B_y_A*B_y_B
    
        
        B_z_A = u_0*I/(4*cp.pi)*(L_x*r_y - L_y*r_x)/((L_x**2 + L_y**2 + L_z**2)*(r_x**2 + r_y**2 + r_z**2) - (r_x*L_x + r_y*L_y + r_z*L_z)**2)
        B_z_B = (L_x**2 + L_y**2 + L_z**2 - r_x*L_x - r_y*L_y - r_z*L_z)/cp.sqrt(r_x**2 + r_y**2 + r_z**2 + 2*(r_x*L_x + r_y*L_y + r_z*L_z) + L_x**2 + L_y**2 + L_z**2) + (r_x*L_x + r_y*L_y + r_z*L_z)/cp.sqrt(r_x**2 + r_y**2 + r_z**2)
        B_z = B_z_A*B_z_B
        
            
        B += cp.array([B_x, B_y, B_z])
        
    return B





# use coordinate transformation in order to set coordinate system to the middle 
X = cp.arange(-0.5*Nx*dx, 0.5*Nx*dx, dx)
Y = cp.arange(-0.5*Ny*dy, 0.5*Ny*dy, dy)
Z = cp.arange(-0.5*Nz*dz, 0.5*Nz*dz, dz)

XX, YY = cp.meshgrid(X, Y)

X, Y, Z = cp.meshgrid(X, Y, Z)


# for a, for now, an arbitrary value as been chosen
B1 = B(X, Y, Z, 0.7*Nz*dz, 0.35*Nz*dz, 10**3, cp.pi/6)
B2 = B(X, Y, Z, -0.7*Nz*dz, 0.35*Nz*dz, -10**3, cp.pi/6)
B3 = Borth(X, Y, Z, 0.7*Nx*dx, 0.35*Nx*dx, 10**3, cp.pi/3)
B4 = Borth(X, Y, Z, -0.7*Nx*dx, 0.35*Nx*dx, -10**3, cp.pi/3)


B = cp.asnumpy(B1 + B2 + B3 + B4)

cp.save('Quad', B)




BU = cp.asnumpy(B1[0] + B2[0] + B3[0] + B4[0])
BV = cp.asnumpy(B1[1] + B2[1] + B3[1] + B4[1])
BZ = cp.asnumpy(B1[2] + B2[2] + B3[2] + B4[2])

BC = cp.asnumpy((BU**2 + BV**2 + BZ**2)**0.5)

fig = plt.figure()
fig.set_size_inches(23, 20)

gs = gridspec.GridSpec(nrows=3, ncols=2, height_ratios=[1, 1, 1])        
ax1 = fig.add_subplot(gs[0, 1])

# set z axis as first axis instead of x
strm = ax1.streamplot(cp.asnumpy(XX), cp.asnumpy(YY),  cp.swapaxes(cp.swapaxes(BU, 0, 2)[30], 0, 1), cp.swapaxes(cp.swapaxes(BV, 0, 2)[30], 0, 1), density=[3, 3], color=cp.swapaxes(BC, 0, 2)[30], linewidth=1, cmap='magma')

fig.colorbar(strm.lines)
ax1.set_clim(vmin=0, vmax=1)
ax1.set_title('Magnetic Field at y = 0')
#plt.savefig('/home/sergey/minerva/python/Quadrupole.png')
