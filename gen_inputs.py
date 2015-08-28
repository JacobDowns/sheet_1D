"""
Generate model inputs for the 1D sheet model
"""

from dolfin import *
from constants import *

# Directory to write model inputs
out_dir = "inputs/"
# Ice sheet length
L = 60e3
# Maximum ice thickness
H_max = 1500.0
# Create the mesh
mesh = IntervalMesh(240, 0.0, L)
V_cg = FunctionSpace(mesh, "CG", 1)


### Create fields

# Bed elevation
B = Function(V_cg)

# Ice thickness
class Thickness(Expression):
    def eval(self,value,x):
        value[0] = sqrt((x[0] + 20.0) * H_max**2 / L)

H = project(Thickness(), V_cg)

# Melt
m = project(Expression("(1.0 + 1.5 * (60000.0 - x[0]) / 60000.0) / 31536000.0"), V_cg)

# Sliding speed
u_b = project(Expression("(10.0 + 240.0 * (60000.0 - x[0]) / 60000.0) / 31536000.0"), V_cg)

# Potential at 0 pressure
phi_m = project(rho_w * g * B, V_cg)

# Overburden pressure
p_i = project(rho_i * g * H, V_cg)

# Potential at overburden pressure
phi_0 = project(phi_m + p_i, V_cg)


### Write fields to a file

File(out_dir + "mesh.xml") << mesh
File(out_dir + "B.xml") << B
File(out_dir + "B.pvd") << B
File(out_dir + "H.xml") << H
File(out_dir + "H.pvd") << H
File(out_dir + "m.xml") << m
File(out_dir + "m.pvd") << m
File(out_dir + "u_b.xml") << u_b
File(out_dir + "u_b.pvd") << u_b
File(out_dir + "phi_m.xml") << phi_m
File(out_dir + "phi_m.pvd") << phi_m
File(out_dir + "p_i.xml") << p_i
File(out_dir + "p_i.pvd") << p_i
File(out_dir + "phi_0.xml") << phi_0
File(out_dir + "phi_0.pvd") << phi_0

"""
plot(B, interactive = True)
plot(H, interactive = True)
plot(m, interactive = True)
plot(u_b, interactive = True)
plot(phi_m, interactive = True)
plot(p_i, interactive = True)
plot(phi_0, interactive = True)
"""
