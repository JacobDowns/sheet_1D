### Create a dictionary of physical constants

physical_constants = {}

# Seconds per day
physical_constants['spd'] = 60.0 * 60.0 * 24.0
# Seconds in a year
physical_constants['spy'] = 60.0 * 60.0 * 24.0 * 365.0                    
# Density of water (kg / m^3)
physical_constants['rho_w'] = 1000.0  
# Density of ice (kg / m^3)
physical_constants['rho_i'] = 910.0
# Gravitational acceleration (m / s^2)
physical_constants['g'] = 9.81 
# Flow rate factor of ice (1 / Pa^3 * s) 
physical_constants['A'] = 5.0e-25
# Average bump height (m)
physical_constants['h_r'] = 0.1
# Typical spacing between bumps (m)
physical_constants['l_r'] = 8.0
# Sheet width under channel (m)
physical_constants['l_c'] = 2.0          
# Sheet conductivity (m^(7/4) / kg^(1/2))
physical_constants['k'] = 1e-2
#k = 5e-3
# Channel conductivity (m^(7/4) / kg^(1/2))
physical_constants['k_c'] = 1e-2
# Specific heat capacity of ice (J / (kg * K))
physical_constants['c_w'] = 4.22e3
# Pressure melting coefficient (J / (kg * K))
physical_constants['c_t'] = 7.5e-8
# Latent heat (J / kg)
physical_constants['L'] = 3.34e5
# Exponents 
physical_constants['alpha'] = 5. / 4.
physical_constants['beta'] = 3. / 2.
physical_constants['delta'] = physical_constants['beta'] - 2.0

