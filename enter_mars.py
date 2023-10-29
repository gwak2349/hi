import numpy as np
import scipy as sp

#M1 is the primary body
#M2 is the secondary body
#R is the distance between the 2 bodies

def determine_sphere_of_influence(M1, M2, R): 
    r_soi = R*(M2/M1)**(2/5)
    return r_soi



def approach_hyperbola(mu,v_infinity,rp):
    eh = 1+(rp*v_infinity**2)/mu
    h = rp*np.sqrt(v_infinity**2+2*mu/rp)
    aiming_radius = h**2/mu*1/np.sqrt(eh**2-1)
    ep = (2*mu-rp*v_infinity**2)/(2*mu+rp*v_infinity**2)
    angle_to_periapsis = np.arccos(1/eh)

    vp1 = np.sqrt(v_infinity**2+2*mu/rp)
    vp2 = np.sqrt(mu*(1+ep)/rp)

    delta_v = vp2-vp1
    return delta_v, aiming_radius, angle_to_periapsis, eh, ep, h



def heliocentric_to_planet(r_sun, v_sun, epsilon):
    
    r1 = [1, 0, 0]
    r2 = [0, np.cos(epsilon), np.sin(epsilon)]
    r3 = [0, -np.sin(epsilon), np.cos(epsilon)]
    matrix = np.array([r1, r2, r3])
    helio_to_planet = np.linalg.inv(matrix)
    
    position = np.array(np.matmul(helio_to_planet,r_sun))

    velocity = np.array(np.matmul(helio_to_planet,v_sun))
    
    return position, velocity

def planet_to_heliocentric(r_mars, v_mars, epsilon):

    r1 = [1, 0, 0]
    r2 = [0, np.cos(epsilon), np.sin(epsilon)]
    r3 = [0, -np.sin(epsilon), np.cos(epsilon)]
    planet_to_helio = np.array([r1, r2, r3])

    ecliptic_position = np.array(np.matmul(planet_to_helio,r_mars))
    ecliptic_velocity = np.array(np.matmul(planet_to_helio,v_mars))

    return ecliptic_position, ecliptic_velocity

def determine_time(h, theta, e, mu):
    F = 2*np.arctanh(np.sqrt((e-1)/(e+1))*np.tan(theta/2))
    M_h = e*np.sinh(F)-F
    t = h**3/mu**2*1/(e**2-1)**(3/2)*M_h

    return M_h, t


def calculate_true_anomaly(e, t_sec, mu, h, initial_guess): #true anomaly for hyperbolic trajectory
    # true_anomaly = []
    t_f = t_sec[-1]
    M_e = mu**2/h**3*(e**2-1)**(3/2)*t_f
    F = initial_guess
    M_h = e*np.sinh(F)-F
    fH = e*np.sinh(F)-F-M_e
    fprime = e*np.cosh(F)-1

    n = 100
    F_0 = M_e
    tol = 1e-8
    for i in range(1,n):
        F_1 = F_0-(fH/fprime)
        F_0 = F_1
        fE = F_0-(e*np.sin(F_0))-M_e
        f_prime = 1-e*np.cos(F_0)
        if abs(F_1-F_0) < tol:
            F = F_1
        else:
            F_1 = F_0


    A = np.tanh(F/2)
    B = np.sqrt((e+1)/(e-1))
    true_anomaly = 2*np.arctan(A/B)
    # true_anomaly.append(ta)

    return true_anomaly




def get_chi(a, h, mu, e, t, r0, theta_0):
    vr0 = mu/h*e*np.sin(theta_0)
    alpha = 1/a
    chi0 = np.sqrt(mu)*abs(alpha)*t
    z = alpha*chi0**2
    
    sz = (np.sinh(np.sqrt(-z))-np.sqrt(-z))/((np.sqrt(-z)))**3
    cz = (np.cosh(-z)-1)/-z
    tol = 1e-8
    ratio = 1
    n = 10000


    def function(chi):
        return r0*vr0/np.sqrt(mu)*chi**2*cz+(1-alpha*r0)*chi**3*sz+r0*chi-np.sqrt(mu)*t
    
    initial_guess = chi0

    root = sp.optimize.fsolve(function, initial_guess, xtol=1e-08)
   

    # fi = r0*vr0/np.sqrt(mu)*chi0**2*cz+(1-alpha*r0)*chi0**3*sz+r0*chi0-np.sqrt(mu)*t
    # fi_prime = r0*vr0/np.sqrt(mu)*chi0*(1-alpha*chi0**2*sz)+(1-alpha*r0)*chi0**2*cz+r0


    
    # for i in range(1,n):
    #     ratio = fi/fi_prime
    #     chi_i = chi0

    #     fi = r0*vr0/np.sqrt(mu)*chi_i**2*cz+(1-alpha*r0)*chi_i**3*sz+r0*chi_i-np.sqrt(mu)*t
    #     fi_prime = r0*vr0/np.sqrt(mu)*chi_i*(1-alpha*chi_i**2*sz)+(1-alpha*r0)*chi_i**2*cz+r0
    #     ratio = fi/fi_prime

        
    #     print(ratio)
    #     if abs(ratio) < tol:
    #         chi = chi_i
    #         chi_list.append(chi)  
    #     else:
    #         chi_i - chi_i-ratio
    # print(chi_list)
    return root[0]
def get_true_anomaly(F0, chi, a, e):
    
    F = F0+chi/np.sqrt(-a)
    theta = 2*np.arctan(np.sqrt((e+1)/(e-1))*np.tanh(F/2))
    return theta


   
def determine_coordinates(r0, v0, delta_theta, mu):
    r = np.linalg.norm(r0)
    v = np.linalg.norm(v0)
    v_r0 = np.dot(r0,v0)/r
    
    h = r*np.sqrt(v**2-v_r0**2)
    
    denominator = 1+(h**2/(mu*r)-1)*np.cos(delta_theta) - h*v_r0/mu*np.sin(delta_theta)
    position = h**2/mu*1/denominator

    f = 1-(mu*position)/h**2*(1-np.cos(delta_theta))
    f_dot = mu/h*(1-np.cos(delta_theta)/np.sin(delta_theta))*(mu/h**2*(1-np.cos(delta_theta))-1/r-1/position)

    g = position*r/h*np.sin(delta_theta)
    g_dot = 1-mu*r/h**2*(1-np.cos(delta_theta))
    
    r_vector = []
    v_vector = []
    for r in range(len(r0)):
        ri = f*r0[r]+g*v0[r]
        r_vector.append(ri)
        vi = f_dot*r0[r]+g_dot*v0[r]
        v_vector.append(vi)


    return r_vector, v_vector

def b_plane_targeting(r0, v0, mu):

    r = np.linalg.norm(r0)
    v = np.linalg.norm(v0)

    h0 = np.cross(r0, v0)/np.linalg.norm(np.cross(r0,v0))
    h = np.linalg.norm(h0)
    e0 = 1/mu*((v**2-mu/r)*r0-(np.dot(r0,v0))*v0)
    e = np.linalg.norm(e0)

    beta = np.arccos(1/e)
    K = np.array([0,0,1])

    k = np.transpose(K)

    S = np.cos(beta)*e0/e+np.sin(beta)*(np.cross(h0,e0)/np.linalg.norm(np.cross(h0,e0)))
    T = np.cross(S,k)/np.linalg.norm(np.cross(S,k))
    R = np.cross(S,T)
    
    a = h**2/mu*1/(e**2-1)
    b  = a*np.sqrt(e**2-1)

    B = b*(np.cross(S,h0))
    BT = np.dot(B,T)
    BR = np.dot(B,R)

    angle = np.arctan(BR/BT)



    # v_infinity, mu, rp


    # B = mu/rp*((1+(v_infinity**2*rp/mu)**2-1)**(0.5)
    
    return B, angle

# https://archive.org/details/explanatorysuppl00pken/page/556/mode/2up


# def rodrigues_rotation_matrix(u, n, alpha):
#     v = (1-np.cos(alpha))*np.dot(u,n)*n+np.cos(alpha)*u-np.sin(alpha)*np.cross(n,u)
#     return v

def calculate_argp(beta, v0):
    vx = v0[0]
    vy = v0[1]
    angle = np.arctan2(vy,vx)
    v0_u = v0/np.linalg.norm(v0)
    x = [1,0,0]
    
    p1 = [np.cos(angle),-np.sin(angle),0]
    p2 = [np.sin(angle),np.cos(angle),0]
    p3 = [0,0,0]
    rotation_matrix = np.array([p1,p2,p3])

    N_vector = np.matmul(x, rotation_matrix)
    
    r1 = [np.cos(beta),np.sin(beta),0]
    r2 = [-np.sin(beta),np.cos(beta),0]
    r3 = [0,0,1]
    rotation_matrix1 = np.array([r1,r2,r3])
    e_vector = np.matmul(v0_u,rotation_matrix1)
    e_unit_vector = e_vector/np.linalg.norm(e_vector)

    N_unit_vector = N_vector/np.linalg.norm(N_vector)
    argp = np.arccos(np.dot(N_unit_vector,e_unit_vector)) #magnitude of vectors being used are 1
    return argp, e_vector
