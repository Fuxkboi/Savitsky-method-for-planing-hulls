import math

LCG = float(input("LCG in [m]: "))
B = float(input("Breadth in [m]: "))
VS = float(input("Ship velocity in [kn]: "))
W = float(input("Ship displacement in [kp]: "))
beta = float(input("Bottom rise angle in [deg]: "))
lamda_correction = float(input("Δλ: "))
friction_coefficient_correction = float(input("ΔCF: "))
BF_correction = float(input("Incorporate Blount - Fox coefficient: "))
accessories_correction = float(input("Accessories correction in [%]: "))


def drag(longitudinal_buoyancy_centre, breadth, ship_velocity, ship_displacement,
         bottom_rise_angle, delta_lamda, delta_cf, bf_coefficient, accessories_coefficient):
    # Gravity acceleration [m/s^2]
    g = 9.81
    # Water density [kp*s^2/m^4]
    rho = 104.61
    # Water kinematic viscosity [m^2/s]
    nu = 1.18831 * 10 ** (-6)
    # Velocity unit transformation ([kn] to [m/s])
    ship_velocity = 0.51444444 * ship_velocity
    # Tolerance used in recursive methods
    recursive_tolerance = 10 ** (-5)
    # Froude number with respect to breadth
    froude_number = ship_velocity / math.sqrt(g * breadth)
    # Average wetted length ratio
    lamda_zero = 4 * longitudinal_buoyancy_centre / (3 * breadth)
    lamda = 0
    flag = True
    while flag:
        lamda = (longitudinal_buoyancy_centre / breadth) / \
                (0.75 - 1 / (5.21 * (froude_number / lamda_zero) ** 2 + 2.39))
        if lamda - lamda_zero < recursive_tolerance:
            flag = False
        else:
            lamda_zero = lamda
    # Prismatic hull lift coefficient
    clb = 2 * ship_displacement / (rho * (ship_velocity ** 2) * (breadth ** 2))
    # Flat plate lift coefficient
    cl0_zero = clb
    cl0 = 0
    flag = True
    while flag:
        cl0 = clb + 0.0065 * bottom_rise_angle * (cl0_zero ** 0.6)
        if cl0 - cl0_zero < recursive_tolerance:
            flag = False
        else:
            cl0_zero = cl0
    # Dynamic trim [deg]
    trim = (cl0 / (0.012 * math.sqrt(lamda) + 0.0055 * (lamda ** 2.5) / froude_number ** 2)) ** (1 / 1.1)
    # Dynamic lift coefficient
    cld = 0.012 * math.sqrt(lamda) * (trim ** 1.1) - 0.0065 * bottom_rise_angle * (
            0.012 * math.sqrt(lamda) * (trim ** 1.1)) ** 0.6
    # Average velocity at hull bottom [m/s]
    average_velocity = ship_velocity * math.sqrt(1 - cld / (lamda * math.cos(trim * math.pi / 180)))
    # Wetted friction length ratio
    lamda_f = lamda + delta_lamda
    # Reynolds number
    reynolds_number = average_velocity * breadth * lamda_f / nu
    # Friction coefficient (ITTC 1957)
    cf = 0.075 / (math.log(reynolds_number, 10) - 2) ** 2 + delta_cf
    # Wetted surface [m^2]
    wetted_surface = lamda_f * (breadth ** 2) / math.cos(bottom_rise_angle * math.pi / 180)
    # Friction component of drag [kp]
    friction_force = 0.5 * rho * (average_velocity ** 2) * wetted_surface * cf
    # Drag force (taking into account given accessories coefficient)
    drag_force = (1 + accessories_coefficient / 100) * \
                 (ship_displacement * math.tan(trim * math.pi / 180) + friction_force / math.cos(trim * math.pi / 180))
    # Blount - Fox coefficient
    displacement_volume = -1
    volumetric_froude_number = -1
    bf = -1
    if bf_coefficient:
        # Displacement volume
        displacement_volume = ship_displacement / (rho * g)
        # Froude number with respect to cube root of displacement volume
        volumetric_froude_number = ship_velocity / math.sqrt(g * displacement_volume ** (1 / 3))
        bf = 0.98 + 2 * ((longitudinal_buoyancy_centre / breadth) ** 1.45) \
            * math.exp(-2 * (volumetric_froude_number - 0.85)) \
            - 3 * (longitudinal_buoyancy_centre / breadth) \
            * math.exp(-3 * (volumetric_froude_number - 0.85))
        drag_force = drag_force * bf
    propulsion_power_demand = 0.00980665 * ship_velocity * drag_force
    return [ship_velocity, froude_number, lamda, clb, cl0, trim,
            cld, average_velocity, lamda_f, reynolds_number, cf,
            wetted_surface, friction_force, displacement_volume, volumetric_froude_number,
            bf, drag_force, propulsion_power_demand]


result = drag(LCG, B, VS, W, beta, lamda_correction, friction_coefficient_correction,
              BF_correction, accessories_correction)
print("Vs " + str(result[0]) + " [m/s]")
print("Fr is " + str(result[1]))
print("λ is " + str(result[2]))
print("CLβ is " + str(result[3]))
print("CL0 is " + str(result[4]))
print("τ is " + str(result[5]) + " [deg]")
print("CLd is " + str(result[6]))
print("Vm " + str(result[7]) + " [m/s]")
print("λF is " + str(result[8]))
print("Re is " + str(result[9]))
print("CF is " + str(result[10]))
print("SF is " + str(result[11]) + " [m^2]")
print("DF is " + str(result[12]) + " [kp]")
print("∇ is " + str(result[13]) + " [m^3]")
print("Fr∇ is " + str(result[14]))
print("M is " + str(result[15]))
print("D is " + str(result[16]) + " [kp]")
print("EHP is " + str(result[17]) + " [kW]")
