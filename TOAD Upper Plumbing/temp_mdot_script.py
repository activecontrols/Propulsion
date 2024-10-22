def pres_mdot_from_prop_mdot(pres_fl, prop_fl, prop_rho, T_tank, P_COPV, prop_mdot):  # Not sure if this will be necessary, but useful just in case??
  # pres_fl: pressurizing fluid (name, string)
  # prop_fl: propellant fluid (name, string)
  # T_tank: pressurizing fluid (K)
  # prop_mdot: propellant mass flow rate (kg/s)
  # returns He_mdot: mass flow rate of Helium needed to get the desired mass
  # flow rate of propellant
  P_COPV = P_COPV * 6894.7572 # psi to Pa
  He_rho = PropsSI('D', 'T', T_tank, 'P', P_COPV, 'Nitrogen')
  He_mdot = prop_mdot * prop_rho / He_rho
  
  return He_mdot

pres_mdot_from_prop_mdot("Helium", "Nitrogen", FUEL_DENSITY, 300, 4500, 10)