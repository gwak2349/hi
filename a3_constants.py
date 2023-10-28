# Constants to be used for AERO2705 Assignments


# -------- STANDARD CONSTANTS --------
G = 6.6726e-11              # Gravitational constant (m^3*kg^(-1)*s^(-2))

Me = 5.974e24               # Mass of Earth (kg)
Re = 6378000                # Radius of Earth (m)
g0 = 9.80665                # Standard sea-level gravity on Earth (m/s^2)
MUe = G*Me                  # Gravitational parameter for earth (m^3*s^(-2))
J2e = 1.08263e-3            # J2 coefficient for Earth
omega_dot = 7.29211510e-5   # Rotation rate of Earth (rad/s)
SOIe = 925e6                # Earth Sphere of Influence (m)
O_Re = 149.6e9              # Earth orbit semimajor axis (m)

Ms = 1.989e30               # Mass of Sun (kg)
MUs = G*Ms                  # Gravitational parameter for Sun (m^3*s^(-2))
Rs = 696e6                  # Radius of Sun (m)

Rm = 3396e3                 # Radius of Mars (m)
Mm = 641.9e21               # Mass of Mars (kg)
MUm = G*Mm                  # Gravitational parameter for Mars (m^3*s^(-2))
J2m = 1.96045e-3            # J2 coefficient for Mars

# Julian day number of current J2000 Julian Epoch (Noon, Jan 1st, 2000)
Julian_epoch = 2451545 

# Define global variables to be used for computations in geocentric, 
# heliocentric frames. 
# The definition of these variables are updated according to the appropriate 
# frame.
M = Me                      # Initialise as Earth mass
MU = MUe                    # Initialise as Earth gravitational parameter
R = Re                      # Initialise as Earth radius
J2 = J2e                    # Initialise J2 perturbation coefficient for Earth


# -------- PROPELLANT CONSTANTS ---------
i_sp_h = 230     # Specific impulse of Monopropellant hydrazine (s)


# -------- TRANSFER CONSTANTS ------------
# Given
r_leo1 = 500e3  # Altitude (m) ---> For Q2
r_leo2 = 600e3   # Altitude (m) --> For Q3
t_atlas = 1549637164021.9214 # ATLAS orbital period




# -------- TYPESETTING CONSTANTS ---------
txtwidth = 6.50127      # Text width of Latex report document