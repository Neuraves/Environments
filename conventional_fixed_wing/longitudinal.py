# longitudinal.py

# Wonhee LEE
# 2024 JUL 05 (FRI)

"""
flight environment for a conventional fixed wing aircraft.
longitudinal dynamics.

- key assumptions
    - incompressible flow
    - troposphere
    - no lateral {forces, moments, winds}
    - no sideslip angle (beta = 0)

- key notations
    - L, D, SF : lift, drag
    - Mm, Lm, Nm : pitching moment, rolling moment, yawing moment
"""

# reference:


import math
import numpy as np


class Aircraft:
    def __init__(self, mass, moi_yy, r_w, r_t, wing, tail, flap, elevator):
        """

        :param mass: total mass
        :param moi_yy: moment of inertia in yy
        :param r_w: r_{wing a.c.} - r_{c.m.}
        :param r_t: r_{tail a.c.} - r_{c.m.}
        """
        self.mass = mass
        self.moi_yy = moi_yy
        self.r_w = r_w
        self.r_t = r_t
        self.wing = wing
        self.tail = tail
        self.flap = flap
        self.elevator = elevator


class Lifting_Surface:
    def __init__(self, C_L0, s, ar, e):
        self.C_L0 = C_L0  # zero-lift lift coefficient
        self.s = s  # span
        self.ar = ar  # aspect ratio
        self.area = s ** 2 / ar
        self.e = e  # Oswald efficiency number
        self.C_l_alpha = 2 * math.pi  # airfoil lift curve slope
        self.C_L_alpha = self.C_l_alpha / (1 + self.C_l_alpha / (math.pi * e * ar))  # wing lift curve slope


class Control_Surface:
    def __init__(self, s, c):
        self.s = s  # span
        self.c = c  # chord
        self.area = s * c
        self.tau = None  # control surface effectiveness


class Air:
    """
    https://www.grc.nasa.gov/www/k-12/airplane/atmosmet.html
    """

    def __init__(self, density=1.225, model="earth_standard_atmosphere", altitude=None):
        self.density = density
        self.model = model
        self.altitude = altitude

    def set_density(self, air_density=None):
        match self.model:
            case "manual":
                self.density = air_density

            case "earth_standard_atmosphere":
                # troposphere
                T = 15.04 - 0.00649 * self.altitude
                p = 101.29 * ((T + 273.1) / 288.08) ** 5.256
                self.density = p / (0.2869 * (T + 273.1))


class Initial_Condition:
    def __init__(self, v=np.zeros((2, 1)), v_w=np.zeros((2, 1)),):
        # in ground-fixed frame: (approximate) inertial frame
        self.v = v  # ground velocity
        self.v_w = v_w  # wind velocity


class Dynamics:
    def __init__(self, aircraft, air, init_cond, eta=1, g=9.81):
        self.aircraft = aircraft
        self.aircraft.flap.tau = control_surface_effectiveness(aircraft.wing, aircraft.flap)
        self.aircraft.elevator.tau = control_surface_effectiveness(aircraft.tail, aircraft.elevator)
        self.air = air
        self.init_cond = init_cond
        self.v = init_cond.v
        self.v_a = self.v - init_cond.v_w
        self.wing_L = None
        self.tail_L = None
        self.eta = eta  # tail dynamic pressure / wing dynamic pressure
        self.g = g  # gravitational acceleration

    def compute_forces_and_moments(self, aircraft, air, v_w, u_f, u_e):
        # air stream
        self.v_a = self.v - v_w
        wing_q = 0.5 * air.density * np.linalg.norm(self.v_a) ** 2
        tail_q = self.eta * wing_q

        # compute change in lift due to control inputs
        wing_delta_L = wing_q * aircraft.wing.area * aircraft.wing.C_L_alpha * aircraft.flap.tau * u_f
        tail_delta_L = tail_q * aircraft.tail.area * aircraft.tail.C_L_alpha * aircraft.elevator.tau * u_e

        # update lift
        """includes lift change due to wind"""
        wing_L = wing_q * aircraft.wing.area * (aircraft.wing.C_L0 + aircraft.wing.C_L_alpha * wing_alpha) \
                 + wing_delta_L
        tail_L = tail_q * aircraft.tail.area * (aircraft.tail.C_L0 + aircraft.tail.C_L_alpha * tail_alpha) \
                 + tail_delta_L
        L = wing_L + tail_L

        # moment
        Mm = -np.cross(aircraft.r_w, wing_L) + -np.cross(aircraft.r_t, tail_L)


def control_surface_effectiveness(lifting_surface, control_surface):
    """
    doi=10.30958/ajte.5-3-2
    """
    return math.sqrt(0.914 * control_surface.area / lifting_surface.area)
