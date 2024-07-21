# test_render_moment_bench.py

# Wonhee LEE
# 2024 JUL 20 (SAT)

"""
test rendering of moment bench environment.
"""

# reference:


import math
import matplotlib.pyplot as plt
import moment_bench


# design airfoil
C_l0 = 0  # symmetric aifoil (flat plate)
c = 0.4  # [m] chord
s = 2  # [m] span
mass = 1  # [kg]
moi_yy = mass * c ** 2 / 12  # [kg*m^2]
l_cm = c / 7
flap = moment_bench.Control_Surface(c=c / 3, s=0.65 * s)

airfoil = moment_bench.Airfoil(C_l0, c, s, mass, moi_yy, l_cm, flap)

# define dynamics
air = moment_bench.Air(density=1.225)

air_speed = 34  # [m/s] wind tunnel air speed
init_cond = moment_bench.Initial_Condition(v_a=air_speed, alpha=math.radians(5), u_f=math.radians(15))

flap_limits = (math.radians(-20), math.radians(20))

dynamics = moment_bench.Dynamics(airfoil, air, init_cond, flap_limits, g=9.81)

# test
dynamics.render()
plt.show()
