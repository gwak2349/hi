import numpy as np

def hyperbolic_departure(mu_sun, mu, R1, R2, rp):

    term1 = np.sqrt(mu_sun/R1)
    term2 = np.sqrt(2*R2/(R1+R2))-1
    v_infinity = term1 * term2
    eh = 1+rp*v_infinity**2/mu
    h = rp*np.sqrt(v_infinity**2+2*mu/rp)
    v1 = np.sqrt(mu/rp)
    v2 = h/rp
    delta_v = v2-v1
    beta = np.arccos(1/eh)
    return delta_v, beta

M_mars = 0.642*1e24
M_earth = 5.97*1e24
M_sun = 1.89*1e30
G = 6.67*1e-11
R1 = 228*1e9
R2 = 149.6*1e9

mars_tilt = 25.19

radius_of_mars = 3389.5*1e3

mu_sun = G*M_sun
mu_earth = G*M_earth
mu_mars = G*M_mars

rp = radius_of_mars+400e3

delta_v, beta = hyperbolic_departure(mu_sun, mu_mars, R1, R2, rp)

print(delta_v, beta*180/np.pi)

#https://www.physicsforums.com/threads/solving-for-time-in-a-hyperbolic-trajectory.1020874/

#velocity value for gmat