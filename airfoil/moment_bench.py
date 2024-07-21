# moment_bench.py

# Wonhee LEE
# 2024 JUL 20 (SAT)

"""
flight environment for a rectangular flat plate airfoil fixed at the center.

- key assumptions
    - incompressible flow
    - no lateral {forces, moments, winds}
    - no sideslip angle (beta = 0)
    - ignore center of mass change due to control surface deflection
    - ignore moment of inertia change due to control surface deflection

- key notations
    - L, D : lift, drag
    - Mm : pitching moment
"""

# reference:


import math
import numpy as np
import matplotlib.pyplot as plt


class Airfoil:
    def __init__(self, C_l0, c, s, mass, moi_yy, l_cm, flap):
        """

        :param mass: total mass
        :param moi_yy: moment of inertia in yy
        :param l_cm: chord-wise distance from the pivot point (center) to the center of mass.
        """
        self.C_l0 = C_l0  # zero-lift lift coefficient
        self.c = c  # chord
        self.s = s  # span
        self.ar = s / c  # aspect ratio
        self.area = self.c * self.s
        self.C_l_alpha = 2 * math.pi  # airfoil lift curve slope

        self.mass = mass
        self.moi_yy = moi_yy
        self.l_cm = l_cm

        self.flap = flap


class Control_Surface:
    def __init__(self, c, s):
        self.c = c  # chord
        self.s = s  # span
        self.ar = s / c  # aspect ratio
        self.area = c * s
        self.tau = None  # control surface effectiveness


class Air:
    def __init__(self, density=1.225):
        self.density = density


class Initial_Condition:
    def __init__(self, v_a=0, alpha=0, u_f=0):
        self.v_a = v_a  # air velocity
        self.alpha = alpha  # angle of attack
        self.u_f = u_f  # flap deflection angle


class Dynamics:
    def __init__(self, airfoil, air, init_cond, flap_limits, g=9.81):
        self.airfoil = airfoil
        self.airfoil.flap.tau = control_surface_effectiveness(airfoil, airfoil.flap)
        self.air = air
        self.init_cond = init_cond
        self.v_a = init_cond.v_a
        self.alpha = init_cond.alpha
        self.u_f = init_cond.u_f
        self.g = g  # gravitational acceleration
        self.W = airfoil.mass * g
        self.flap_limits = flap_limits  # flap deflection limits (lower bound, upper bound)

        self.L = 0
        self.Mm = 0

    def compute_forces_and_moments(self, airfoil, air, u_f):
        # air stream
        q = 0.5 * air.density * np.linalg.norm(self.v_a) ** 2

        # compute change in lift due to control inputs
        u_f = self.limit_control_input(u_f, self.flap_limits)
        delta_L = q * airfoil.area * airfoil.C_l_alpha * airfoil.flap.tau * u_f
        # print(f"u_f: {math.degrees(u_f):.4F} [deg], delta_L: {delta_L:.4F} [N]")

        # update lift
        """includes lift change due to wind"""
        L = q * airfoil.area * (airfoil.C_l0 + airfoil.C_l_alpha * self.alpha) + delta_L
        # print(f"L: {L:.4F} [N], W: {self.W:.4F} [N]")

        # moment
        arm_ac = airfoil.c / 4 * np.cos(self.alpha)
        arm_cm = airfoil.l_cm * np.cos(self.alpha)
        Mm = arm_ac * L - arm_cm * self.W

        self.L = L
        self.Mm = Mm
        return L, Mm

    def limit_control_input(self, u, limits):
        if u < limits[0]:
            u = limits[0]
        elif u > limits[1]:
            u = limits[1]

        return u

    def reset(self):
        self.v_a = self.init_cond.v_a
        self.alpha = self.init_cond.alpha
        self.u_f = self.init_cond.u_f
        self.L = 0
        self.Mm = 0


def control_surface_effectiveness(airfoil, control_surface):
    """
    doi=10.30958/ajte.5-3-2
    """
    return math.sqrt(0.914 * control_surface.area / airfoil.area)


def R(alpha):
    return np.array([[np.cos(alpha), -np.sin(alpha)],
                     [np.sin(alpha), np.cos(alpha)]])
