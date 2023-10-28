# Contains orbital transform functions for AERO2705 Assignments

import numpy as np
import scipy as sp
import a3_constants as c

def masked_ra(ra):
    # Compute masked array of values for 
    # plotting discontinuities

    # Compute difference between successive values
    ra_diffs = np.append(np.diff(ra), 0)

    # Find differences that are large
    dis_idx = np.abs(ra_diffs) > 1.75*np.pi

    # Mask array at discontinuities
    m_ra = np.ma.masked_array(ra,dis_idx)

    return m_ra

# Compute specific angular momentum given semimajor axis and eccentricity
def h(a, e):
    return (a*c.MU*(1 - e**2))**(1/2)

# Compute semimajor axis given specific angular momentum and eccentricity
def a(h, e):
    return h**2/(c.MU*(1 - e**2))

# Semimajor axis given perigee, apogee radii
def a(p,a):
    return (p + a)/2   

# Eccentricity given perigee, apogee radii
def e(p,a):
    return (a - p)/(a + p) 

def set_geo_const():
    # Set constants to their geocentric values
    c.R = c.Re
    c.MU = c.MUe
    c.M = c.Me
    c.J2 = c.J2e
    return

def set_helio_const():
    # Set constants to heliocentric values
    c.R = c.Rs
    c.MU = c.MUs
    c.M = c.Ms
    return

def set_mars_const():
    # Set constants to mars-centered values
    c.R = c.Rm
    c.MU = c.MUm
    c.M = c.Mm
    c.J2 = c.J2m
    return

# -------- ORBITAL DATA EXTRACTION --------------------------------------

def lead_zero(coeff,exp):
    #Adds a leading decimal for given coefficient & power of 10
    return (
        int(coeff)*10**(int(exp) - (len(str(coeff)) - str(coeff).count(" ")))
        )

def date_time(jul):
    #Returns day, hour, minute, sec from given julian date fraction
    day = np.floor(jul)
    hrd = (jul - day)*24            # Hour decimal
    hr  = np.floor(hrd)
    min = np.floor((hrd - hr)*60)
    sec = (((hrd - hr)*60)-min)*60
    return [day, hr, min, sec]

def tle_extract():
    #Takes TLE as input and returns orbital elements and satellite information
    # in SI units
    # Assumes TLE file has 3 lines and is in correct format

    f = open("EASA_XMM_tle.txt")
    tle_lines = f.readlines()
    f.close()

    #Extract TLE lines (including newline (\n) characters)
    name    = tle_lines[0][:]
    line1   = tle_lines[1][:]
    line2   = tle_lines[2][:]

    #Extract elements from TLE line 1
    cat1    = int(line1[2:7])
    clssf   = line1[7]
    intl_y  = int(line1[9:11])
    intl_n  = int(line1[11:14])
    intl_pc = line1[14:17]
    epoch   = np.append(                # (year, (day, hour, minute, second))
        [int(line1[18:20])], 
        date_time(float(line1[20:32]))
        )   
    fdmm    = float(line1[33:43]) 
    sdmm    = lead_zero(                # Apply leading decimal
        line1[44:50],
        line1[50:52]
          )
    bstar   = lead_zero(                # Apply leading decimal
        line1[53:59], 
        line1[59:61]
        )
    eph     = int(line1[62])
    set_n   = int(line1[64:68])
    chk1    = int(line1[68])

    # Extract elements from TLE line 2
    cat2    = int(line2[2:7])
    i       = np.radians(float(line2[8:16]))            # (rad)
    raan    = np.radians(float(line2[17:25]))           # (rad)
    e       = lead_zero(line2[26:33],0)
    argp    = np.radians(float(line2[34:42]))           # (rad)
    ma      = np.radians(float(line2[43:51]))           # (rad)
    mm      = (2*np.pi*float(line2[52:63]))/(24*60**2)  # (rad/s)
    t       = 2*np.pi/mm                                # Orbital period (s)
    rev_n   = int(line2[63:68])
    chk2    = int(line2[68])

    # Kepler's 3rd law & mean motion used to calculate semimajor axis.
    a = ((t**2*c.MU)/(4*np.pi**2))**(1/3)

    # Define dictionaries to be returned.
    elements = {
        "epoch": epoch,
        "fdmm": fdmm,
        "sdmm": sdmm,
        "bstar": bstar,
        "i": i,
        "raan": raan,
        "e": e,
        "argp": argp,
        "ma": ma,
        "mm": mm,
        "a": a
        }
    sat_info = {
        "name": name,
        "cat_1": cat1,
        "clssf": clssf,
        "intl_y": intl_y,
        "intl_n": intl_n,
        "intl_pc": intl_pc,
        "eph": eph,
        "set_n": set_n,
        "chk_1": chk1,
        "cat_2": cat2,
        "rev_n": rev_n,
        "chk2": chk2
        }
    
    return elements,sat_info

def txt_extract():

    f = open("comet_atlas_data.txt")
    data_lines = f.readlines()
    f.close()

    jul_frac = (
        float(data_lines[3]) 
        - 2459945.5
        )                     # Julian Day fraction (ordinal date) at epoch
    epoch = np.append(
        [int(data_lines[5])],
        date_time(jul_frac)
        )                     # [year, day, hour, minute, second]

    # Define dictionaries to be returned.
    elements = {
        "epoch": epoch,
        "fdmm": 0,
        "sdmm": 0,
        "bstar": 0,
        "perihelion": (int(data_lines[9])*1e3),
        "i": np.radians(float(data_lines[21])),
        "raan": np.radians(float(data_lines[23])),
        "e": float(data_lines[7]),
        "argp": np.radians(float(data_lines[25])),
        "ma": np.radians(float(data_lines[27])),
        "mm": (np.radians(float(data_lines[29])))/(24*60**2),
        "a": (float(data_lines[13])*1e3)
        }
    sat_info = {
        "name": data_lines[1]
        }
    
    return elements,sat_info

def compute_jd(t):
    # Compute Julian Date (JD) given (yy, DOY, hh, mm, ss)

    # Extract data from input
    y = t[0]    # Year - 2000
    d = t[1]    # Day number
    h = t[2]    # Hour
    m = t[3]    # Minute
    s = t[4]    # Second
    a = 0       # 1 if leap year, otherwise 0

    # Compute Julian Day
    julian_day = (
        2451544.5 
        + 365*y
        + np.floor(0.25*y)
        - np.floor(0.01*y)
        + np.floor(0.0025*y)
        + a
        + d
        + (1/24)*h
        + (1/1440)*m
        + (1/86400)*s
        )

    return julian_day


# -------- ORBITAL ELEMENTS --------------------------------------

def compute_anomalies(t_sec, ma, e, mm, hyp):
    # Compute true, eccentric and mean anomalies from mean anomaly, 
    # eccentricity and mean motion

    # Propagate mean anomaly using mean motion (vector over time)
    ma_vec = mm*np.copy(t_sec) + ma
    
    # Define Kepler equation as a local function (for optimisation)
    # Uses Newton's method to find roots of Kepler equation
    def kepler_equation(E):
        
        # Set desired absolute tolerance
        tol = 1e-12

        # Define kepler equation to be called during iteration
        f = lambda E: E - e*np.sin(E) - np.copy(ma_vec)

        # Iteration step (Newton-Raphson method)
        while (np.max(np.absolute(f(E))) > tol):
            E -= f(E)/(1 - e*np.cos(E))
        return E
    
    def hyp_kepler_equation(H):

        # Set desired absolute tolerance
        tol = 1e-12

        # Define kepler equation to be called during iteration
        f = lambda H: -e*np.sinh(H) + H + np.copy(ma_vec)

        # Iteration step (Newton-Raphson method)
        while (np.max(np.absolute(f(H))) > tol):
            H += f(H)/(e*np.cosh(H) - 1)
        return H

    # Define our initial guess for the eccentric anomaly
    init_guess = np.copy(ma_vec)

    if hyp == True:

        # Calc hyperbolic anomaly
        ea = hyp_kepler_equation(init_guess)

        # Calc true anomaly
        ta = 2*np.arctan2(np.tanh(np.copy(ea)/2)*(e + 1)**(1/2),(e - 1)**(1/2))
    else:
        # Calc eccentric anomaly
        ea = kepler_equation(init_guess)

        # Calc true anomaly
        ta = 2*np.arctan2(np.tan(np.copy(ea)/2)*(1 + e)**(1/2),(1 - e)**(1/2))

    # Wrap anomalies to -pi:pi
    ea = np.arctan2(np.sin(ea),np.cos(ta))
    ma_vec = np.arctan2(np.sin(ma_vec),np.cos(ma_vec))

    return ta, ea, ma_vec


def compute_orbital_velocity(a, e, ta):
    # Computes radial and tangential components of velocity given semimajor 
    # axis, eccentricity, true anomaly and gravitational parameter

    v_r = c.MU/h(a, e)*e*np.sin(np.copy(ta)) # Radial velocity
    v_n = c.MU/h(a, e)*(1 + e*np.cos(np.copy(ta))) # Normal/Tangential velocity

    return v_r, v_n

def compute_r(a, e, ta):
    # Computes magnitude of position given semimajor axis, eccentricity and 
    # true anomaly
    
    r_vec = a*(1 - e**2)/(1 + e*np.cos(np.copy(ta)))

    return r_vec

def compute_flightpath_angle(e, theta):
    # Compute flightpath angle relative to local horizon given eccentricity 
    # and true anomaly

    gamma = np.arctan2(e*np.sin(theta),1 + e*np.cos(theta))

    return gamma

def compute_ma(period, t):
    # Compute mean anomaly given an orbital period and a time t
    return 2*np.pi*t/period

def compute_theta(ep):
    # Compute Greenwich sidereal time over duration of simulation from 
    # julian date at epoch

    # Compute Julian day at epoch   
    jd = compute_jd(ep)

    # Get julian day at 0h UT
    j0 = np.floor(jd) + 0.5

    # Compute Julian centuries fraction since J2000 epoch
    t0 = (j0 - 2451545)/36525
    
    # Compute Greenwich sidereal time at Oh UT (rad) and 
    # reduce into range (0, 2pi)
    theta_g0 = np.radians(
        100.4606184
        + 36000.77004*t0
        + 0.000387933*t0**2
        - 2.583e-8*t0**3
        ) % 2*np.pi
    
    # Compute Greenwich sidereal time at epoch
    theta_g = (theta_g0 + np.radians(360.98564724)*(jd - j0)) % 2*np.pi

    return theta_g

# ----------- HYPERBOLIC TRAJECTORIES ---------------------------------------

def compute_e_hyp(rp,v_inf):
    # Compute eccentricity of hyperbolic orbit given periapsis radius 
    # and hyperbolic excess speed
    return 1 + (rp*v_inf**2)/c.MU

def compute_h_hyp(rp,v_inf):
    # Compute specific angular momentum of hyperbolic orbit
    return  rp*np.sqrt(v_inf**2 + (2*c.MU)/rp)

def compute_a_hyp(v_inf):
    # Compute hyperbolic semimajor axis
    return - c.MU/v_inf**2

def compute_mm_hyp(a_hyp):
    # Compute hyperbolic mean motion
    return np.sqrt(c.MU/(-a_hyp)**3)

def compute_circ_orbit_v(rp):
    # Velocity of circular parking orbit given periapse radius
    return np.sqrt(c.MU/rp)

def compute_escape_deltav(v_hyp_perigee,vc):
    # Compute escape delta V required for hyperbolic escape trajectory
    return v_hyp_perigee - vc

def compute_hyp_conic_parameter(a_hyp,e_hyp):
    # Compute conic parameter given semimajor axis and eccentricity
    return a_hyp*(1 - e_hyp**2)

def compute_ta_hyp(a,e,r):
    # Compute true anomaly for given radius, eccentricity and semimajor axis
    return np.arccos((a*(1-e**2) - r)/(e*r))

def compute_ha(e,ta):
    # Compute hyperbolic anomaly given eccentricity and true anomaly
    return 2*np.arctan2(np.sqrt(e - 1)*np.tan(ta/2),np.sqrt(e + 1))

def compute_hyp_ma(e,ha):
    # Compute hyperbolic mean anomaly given eccentricity and hyperbolic anomaly
    return e*np.sinh(ha) - ha
    

# -------- ORBITAL PERTURBATIONS --------------------------------------

def compute_j2(j2_coeff, R, a, e, i):
    # Compute J2 perturbations for a given orbit. 

    # Compute node-line regression rate
    raan_dot = -(
        3*np.sqrt(c.MU)*j2_coeff*R**2
        )/(
            2*(1 - e**2)**2*a**(7/2)
            )*np.cos(i)

    # Compute rate of perigee advance
    argp_dot = -(
        3*np.sqrt(c.MU)*j2_coeff*R**2
        )/(
            2*(1 - e**2)**2*a**(7/2)
            )*(
                (5/2)*(np.sin(i))**2 - 2
                )

    return raan_dot, argp_dot


# -------- COORDINATE TRANSFORMATIONS --------------------------------------

def elements_to_perifocal(ta, a, e, mu, h):
    # Calc perifocal distance in m
    r = compute_r(a, e, ta)

    # Compute perifocal coordinates
    p = np.copy(r)*np.cos(np.copy(ta))
    q = np.copy(r)*np.sin(np.copy(ta))
    w = np.zeros_like(p)

    # Compute perifocal velocities
    dp = - mu/h*np.sin(np.copy(ta))
    dq = mu/h*(e + np.cos(np.copy(ta)))
    dw = np.zeros_like(dp)

    return p, q, w, dp, dq, dw

def perifocal_to_eci(p, q, w, dp, dq, dw, i, raan, argp):
    # Transform coordinates to ECI frame

    # Initialise arrays for transformed position and velocity components
    x = np.zeros_like(p)
    y = np.zeros_like(p)
    z = np.zeros_like(p)

    dx = np.zeros_like(p)
    dy = np.zeros_like(p)
    dz = np.zeros_like(p)

    # Create position vector
    r_perifocal = np.array([p, q, w])

    # Compute transformation of position vector
    r_eci = np.einsum(
        'ijk,jk->ik',
        perifocal_to_eci_matrix(i, raan, argp),
        r_perifocal
        )

    # Extract transformed position components
    x = r_eci[0]
    y = r_eci[1]
    z = r_eci[2]

    # Create velocity vector
    rdot_perifocal = np.array([dp, dq, dw])

    # Compute transformation of velcity vector
    rdot_eci = np.einsum(
        'ijk,jk->ik',
        perifocal_to_eci_matrix(i, raan, argp),
        rdot_perifocal
        )

    # Extract transformed velocity components
    dx = rdot_eci[0]
    dy = rdot_eci[1]
    dz = rdot_eci[2]

    return x, y, z, dx, dy, dz

def eci_to_ecef(x, y, z, dx, dy, dz, theta):
    # Transform coordinates to ECEF frame

    # Initialise arrays for transformed position and velocity components
    x_dash = np.zeros_like(x)
    y_dash = np.zeros_like(x)
    z_dash = np.zeros_like(x)
    dx_dash = np.zeros_like(x)
    dy_dash = np.zeros_like(x)
    dz_dash = np.zeros_like(x)

    # Create vector of ECI position components
    r_eci = np.array([x, y, z])

    # Apply transformation matrix to position vector
    r_ecef = np.einsum(
        'ijk,jk->ik',
        eci_to_ecef_matrix(theta),
        r_eci
        )
    
    # Extract individual ECEF position components
    x_dash = r_ecef[0]
    y_dash = r_ecef[1]
    z_dash = r_ecef[2]

    # Create vector of ECI velocity components
    v_eci = np.array([dx, dy, dz])

    # Apply transformation matrix to velocity vector
    v_ecef = np.einsum(
        'ijk,jk->ik',
        eci_to_ecef_matrix(theta),
        v_eci
        )
    
    # Extract individual ECEF velocity components
    dx_dash = v_ecef[0]
    dy_dash = v_ecef[1]
    dz_dash = v_ecef[2]

    return x_dash, y_dash, z_dash, dx_dash, dy_dash, dz_dash

def ra_and_dec(x_dash, y_dash, z_dash):
    # Compute right ascension and declination from coordinates

    # Compute magnitude of r
    r = np.sqrt(x_dash**2 + y_dash**2 + z_dash**2)

    # Compute direction cosines
    l = x_dash/r
    m = y_dash/r
    n = z_dash/r

    # Compute declination
    dec = np.arcsin(n)

    # Define function to compute right ascension
    def compute_ra(m, l, dec):
        # Compute right ascension
        if m > 0:
            ra = np.arccos(l/np.cos(dec))
        elif m <= 0:
            ra = 2*np.pi - np.arccos(l/np.cos(dec))
        return ra
    
    vcompute_ra = np.vectorize(compute_ra)

    ra = vcompute_ra(m, l, dec)

    return ra, dec


def heliocentric_to_planet(vec, epsilon):
    # Coordinate transform from heliocentric to planet-centered inertial
    # Epsilon is the angle between the helio NP and planet NP

    r1 = [1, 0, 0]
    r2 = [0, np.cos(epsilon), np.sin(epsilon)]
    r3 = [0, -np.sin(epsilon), np.cos(epsilon)]
    matrix = np.array([r1, r2, r3])
    helio_to_planet = np.linalg.inv(matrix)

    return np.array(np.matmul(helio_to_planet,vec))

def planet_to_heliocentric(vec, epsilon):
    # Coordinate transform from planet-centered inertial to heliocentric
    # Epsilon is the angle between the helio NP and planet NP

    r1 = [1, 0, 0]
    r2 = [0, np.cos(epsilon), np.sin(epsilon)]
    r3 = [0, -np.sin(epsilon), np.cos(epsilon)]
    planet_to_helio = np.array([r1, r2, r3])

    return np.array(np.matmul(planet_to_helio,vec))



# -------- COORDINATE TRANSFORMATION MATRICES --------------------------------

def perifocal_to_eci_matrix(i, raan, argp):
    # Calculate transformation matrix from perifocal to ECI frame

    p_to_e = [
        [-np.sin(raan)*np.cos(i)*np.sin(argp) + np.cos(raan)*np.cos(argp),
         -np.sin(raan)*np.cos(i)*np.cos(argp) - np.cos(raan)*np.sin(argp),
         np.sin(raan)*np.sin(i)],
        [np.cos(raan)*np.cos(i)*np.sin(argp) + np.sin(raan)*np.cos(argp),
         np.cos(raan)*np.cos(i)*np.cos(argp) - np.sin(raan)*np.sin(argp),
         -np.cos(raan)*np.sin(i)],
        [np.sin(i)*np.sin(argp),
         np.sin(i)*np.cos(argp),
         np.cos(i)]
        ]

    return p_to_e

def eci_to_ecef_matrix(theta):
    # Calculate transformation matrix from ECI to ECEF frame

    e_to_ecef = np.array([
        [np.cos(theta), np.sin(theta), np.zeros_like(theta)],
        [-np.sin(theta), np.cos(theta), np.zeros_like(theta)],
        [np.zeros_like(theta), np.zeros_like(theta), np.ones_like(theta)]
    ])

    return e_to_ecef

# STATE VECTOR/KEPLERIAN ELEMENTS CONVERSIONS

def state_vecs_to_elements(ro, vo, mu): #ro and vo are initial position and velocity vectors respectively, 

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
    
    e_vector = 1/mu*((voi**2-mu/roi)*ro-roi*v_r*vo)
    e_x = e_vector[0]
    e_y = e_vector[1]
    e_z = e_vector[2]
    e = np.linalg.norm(e_vector)
    if e_z < 0:
        argp = 2*np.pi - np.arccos(np.dot(N_vector, e_vector)/(N*e))
    else:
        argp = np.arccos(np.dot(N_vector, e_vector)/(N*e))
    
    if v_r < 0:
        theta = 2*np.pi - np.arccos(np.dot(N_vector, e_vector)/(N*e))
    else:
        theta = np.arccos(np.dot(N_vector, e_vector)/(N*e))
    
    return h, i, raan, e, argp, theta 

def elements_to_state_vecs(a,e,i,argp,raan,ta):
    # Compute state vectors for given true anomaly and keplerian elements

    p = a*(1 - e**2)        # Conic parameter
    r = compute_r(a,e,ta)   # Compute radius
    h_orbit = h(a,e)        # Angular momentum

    # Position state vector
    x = r*(np.cos(raan)*np.cos(argp + ta) - np.sin(raan)*np.sin(argp + ta)*np.cos(i))
    y = r*(np.sin(raan)*np.cos(argp + ta) + np.cos(raan)*np.sin(argp + ta)*np.cos(i))
    z = r*(np.sin(i)*np.sin(argp + ta))

    r_vec = np.array([x,y,z])

    # Velocity state vector
    vx = x*h_orbit*e/(r*p)*np.sin(ta) - h_orbit/r*(np.cos(raan)*np.sin(argp + ta) + np.sin(raan)*np.cos(argp + ta)*np.cos(i))
    vy = y*h_orbit*e/(r*p)*np.sin(ta) - h_orbit/r*(np.sin(raan)*np.sin(argp + ta) - np.cos(raan)*np.cos(argp + ta)*np.cos(i))
    vz = z*h_orbit*e/(r*p)*np.sin(ta) + h_orbit/r*np.sin(i)*np.cos(argp + ta)

    r_dot_vec = np.array([vx,vy,vz])

    state_vecs = [r_vec,r_dot_vec]

    return state_vecs




