import numpy as np
import pylab as plt

hc = 197.327 # Conversion factor in MeV fm (hut * c)
G = hc * 6.67259e-45 # Gravitational constant
Ms = 1.1157467e60
density = 1665.3 # Central density (density at r = 0)
M0 = (4*3.14159265*(G**3)*density)**(-0.5)
R0 = G*M0
mn = 938.926 # Mass of neutron in MeV c^-2

def neutron():          
    n = 1
    err = 1
    tol = 1e-15
    count = 0
    # Newton method
    while err > tol : 
        count += 1
        fn = n*mn + 236*n**(2.54) - density
        dfn = mn + 236*2.54*n**(1.54)
        temp = n - fn/dfn
        err = np.abs(n-temp)
        n = temp
    print ("Newton Converged after ", count, "iterations")
    return n
    
    
def density_function(p):
    n = (p*density/363.44)**(1./2.54)
    return (236. * n**2.54 + n *mn)/density 
    

def pressure_radius(r,m,p,flag):
    if flag == 0: # Classical model
        y = -m*density_function(p)/(r**2 + 1e-20)
    else: # Relativistic model
        rh = density_function(p)                            
        y = -(p+rh)*(p*r**3 + m)/(r**2 - 2*m*r + 1e-20)
    return y

def mass_radius(r,m,p):
    return density_function(p)*r**2
    

def Euler(r,m,p,h,flag):
    y = np.zeros(2)
    y[0] = m + mass_radius(r,m,p)*h
    y[1] = p + pressure_radius(r,m,p,flag)*h
    return y
    
def RK4(r,m,p,h,flag):
    y = np.zeros(2)
    a1 = mass_radius(r,m,p)
    b1 = pressure_radius(r,m,p,flag)
    
    a2 = mass_radius(r+0.5*h,m+0.5*a1*h,p+0.5*b1*h)
    b2 = pressure_radius(r+0.5*h,m+0.5*a1*h,p+0.5*b1*h,flag)
    
    a3 = mass_radius(r+0.5*h,m+0.5*a2*h,p+0.5*b2*h)
    b3 = pressure_radius(r+0.5*h,m+0.5*a2*h,p+0.5*b2*h,flag)
    
    a4 = mass_radius(r+h,m+h*a3,p+h*b3)    
    b4 = pressure_radius(r+h,m+h*a3,p+h*b3,flag)    
    
    y[0] = m + h*(a1 + 2.*a2 + 2.*a3 + a4)/6.
    y[1] = p + h*(b1 + 2.*b2 + 2.*b3 + b4)/6.
    return y
    

def mplot(fign,x,y,xl,yl,clr,lbl):
    plt.figure(fign)
    plt.xlabel(xl)    
    plt.ylabel(yl)
    plt.title("Relativistic model")
    return plt.plot(x,y,clr, linewidth =2.0,label = lbl)

plt.figure(1)
plt.clf
plt.figure(2)
plt.clf


N = 4000
r = np.linspace(0,15,N)
h = r[1]-r[0]

m = np.zeros(N)
p = np.zeros(N)
rh = np.zeros(N)
ni = neutron()

r[0]   = 0
m[0] = 0
p[0] = 363.44 * (ni**2.54)/density
rh[0] = 1



flag_set = [0,1]
print ("Initial number density, ni =", ni)
print ("Initial Pressure, P[0] =", p[0])

print ("Simulation range, R = 0 to", r[-1]*R0*1e-18, "km")
tol = 9e-5


# Looping over the two models (k = 0: classical, k = 1: relativistic)
for k in range(0,2):
    flag = flag_set[k]
    for i in range(0,N-1):
        [m[i+1], p[i+1]] = RK4(r[i],m[i],p[i],h,flag)
        rh[i+1] = density_function(r[i])
        if p[i+1] < tol:
            rf = r[i]
            mf = m[i]
            break

    print ("Running RK4")
    if i == N-2:
        print ("Program didn't converge to P = 0, extend the maximum value of r")
    else:
        print ("P <", tol, "found after", i, "runs")
    
    m = m[0:i+2]
    p = p[0:i+2]
    rh = rh[0:i+2]
    r = r[0:i+2]
    
    if flag == 0:
        lbl = "Classical Model"
        clr = "green"
    else:
        lbl = "Relativistic Model"
        clr = "blue"

    print ("----------------------------------------------")
    print (lbl, "Results\n")
    print ("----------------------------------------------")
    print ("Initial density, density =", density, "MeV/fm3")
    print ("Total mass = ", mf*M0/Ms, "times Solar mass")
    print ("Radius of the Neutron star =", rf*R0*1e-18, "km")


    mplot(1,r*R0*1e-18,p*density,'radius (km)','pressure',clr,lbl)
    mplot(2,r*R0*1e-18,m*M0/Ms,'radius (km)','mass (solar mass)',clr,lbl)
    
plt.figure(1)
q = plt.legend(loc = 0)
q.draw_frame(False)
plt.figure(2)
q=plt.legend(loc = 0)
q.draw_frame(False)
