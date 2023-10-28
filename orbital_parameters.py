import numpy as np

def calculate_orbital_parameters(ro, vo, mu): #ro and vo are initial position and velocity vectors respectively, 

    roi = np.linalg.norm(ro)
    voi = np.linalg.norm(vo)
    v_r = np.dot(ro,vo/roi)
    
    
    K = np.array([0,0,1])   
    h_vector = np.cross(ro, vo)
    h_x = h_vector[0]
    h_y = h_vector[1]
    h_z = h_vector[2]
    h = np.linalg.norm(h_vector)

    i = np.arccos(h_z/h)


    N_vector = np.cross(K, h_vector)

    N_x = N_vector[0]
    N_y = N_vector[1]
    N = np.linalg.norm(N_vector)

    if N_y < 0:
        raan = 2*np.pi-np.arccos(N_x/N)
    else:
        raan = np.arccos(N_x/N)
    e_vector = []
    for i in range(len(ro)):
        t1 = voi**2-mu/roi
        t2 = roi*v_r
        evector = 1/mu*(t1*ro[i]-t2*vo[i])
        e_vector.append(evector)
    e_x = e_vector[0]
    e_y = e_vector[1]
    e_z = e_vector[2]
    e = np.linalg.norm(e_vector)
    if e_z < 0:
        argp = 2*np.pi - np.arccos(np.dot(N_vector, e_vector)/(N*e))
    else:
        argp = np.arccos(np.dot(N_vector, e_vector)/(N*e))
    
    if v_r < 0:
        theta = 2*np.pi - np.arccos(np.dot(e_vector, ro)/(roi*e))
    else:
        theta = np.arccos(np.dot(e_vector, ro)/(roi*e))
    
    return h, i, raan, e, argp, theta 




