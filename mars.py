import matplotlib.pyplot as plt
import numpy as np
import enter_mars as em
import orbital_parameters as op
import a3_orbitalTransforms as aot

M_mars = 0.642*1e24
M_earth = 5.97*1e24
M_sun = 1.89*1e30
G = 6.67*1e-11
R1 = 149.6*1e9
R2 = 228*1e9
mars_tilt = 25.19*np.pi/180

radius_of_mars = 3389.5e3

mu_sun = G*M_sun
mu_earth = G*M_earth
mu_mars = G*M_mars

r_soi = em.determine_sphere_of_influence(M_sun, M_mars, R2)


v_infinity = 2526.442

v_mars = [-21263.13376563, 10506.3778946, 5392.26564194]
unit_vector = v_mars/np.linalg.norm(v_mars) 					



r_sun = r_soi*unit_vector
v_sun = v_infinity*unit_vector


rp = 400e3+radius_of_mars
# e = 0.4
# ra = rp*(1+e)/(1-e)

# print(ra-radius_of_mars)



# v_infinity = 3918.6075818730897
delta_v, aiming_radius, angle_to_periapsis, eh, ep, h_hyperbola = em.approach_hyperbola(mu_mars, v_infinity, rp)
print(eh)
# print(ep)
print(delta_v)
# print(h_hyperbola)

# a = h**2/mu_mars*1/(1-ep**2)

r_soi  = em.determine_sphere_of_influence(M_sun, M_mars, R2)

vp = np.sqrt(mu_mars/rp)

h_orbit = rp*vp
ra = h_orbit**2/mu_mars*1/(1-ep)
a1 = h_orbit**2/mu_mars*1/(1-ep**2)

# print(ra, a1)

# distance_to_soi = R2-r_soi


theta_infinity = np.arccos(-1/eh)

# print(theta_infinity*180/np.pi)

epsilon = 23.44*np.pi/180

r0_mars, v0 = em.heliocentric_to_planet(r_sun, v_sun, epsilon)
random_vector = [0,0,1]

r = np.cross(r0_mars, random_vector)
r0 = r0_mars + aiming_radius*r


# h1, i1, raan1, e1, argp1, theta1 = op.calculate_orbital_parameters(r0, v0, mu_mars) 

# if h1>h_hyperbola:
#     print("hi")

# print(h1, h_hyperbola)


# ma_initial, t_final = em.determine_time(h_hyperbola, theta_infinity, eh, mu_mars)



# a_hyperbola = aot.compute_a_hyp(v_infinity)
# F0 = aot.compute_ha(eh,theta_infinity)

# mm = aot.compute_mm_hyp(a_hyperbola)
# ma = aot.compute_hyp_ma(eh,F0)
# t_flight = ma/mm
# t_sec = np.linspace(0,t_flight, 100)

# # print(t_flight, t_final)

# # true_anomaly = em.calculate_true_anomaly(eh, t_sec, mu_mars, h_hyperbola, ma)


# root = em.get_chi(a_hyperbola, h_hyperbola, mu_mars, eh, t_flight, r_soi, theta_infinity)
# true_anomaly = em.get_true_anomaly(F0, root, a_hyperbola, eh,)

# print(true_anomaly)
# r_vector, v_vector = em.determine_coordinates(r0, v0, theta_infinity-true_anomaly, mu_mars)

# perigee = np.linalg.norm(r_vector)
# # h, i, raan, e, argp, theta = op.calculate_orbital_parameters(r_vector, v_vector, mu_mars)


# # print(true_anomaly*180/np.pi)
# print(perigee, rp)
# if perigee == rp:
#     print("nay")

# 6902877174341330.0

# 6902874987106839.0



# # print(r0_mars)
# # print(v0_mars)

# # r_magnitude = np.linalg.norm(r0_mars)
# # distance = R2-r_soi
# # distance1 = R2 - r_magnitude
# # print(distance1, distance)
# # if distance1 < distance:
# #     print("No")

# # r0 = [r_soi*np.cos(declination)*np.cos(right_ascension),r_soi*np.cos(declination)*np.sin(right_ascension), r_soi*np.sin(declination)]


# # rp1, vp = em.determine_coordinates(r0, v0, theta_infinity, mu_mars)



# # r_magnitude = np.linalg.norm(rp1)

# # print(r_magnitude, rp)
# # if r_magnitude == rp:
# #     print("Hi")

# https://www.youtube.com/watch?v=Fh3nMi87cB8 - rotation matrix