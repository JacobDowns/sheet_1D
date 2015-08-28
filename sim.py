"""
1D Sheet model sim
"""

from dolfin import *
from sheet_model import *

# Model input directory
in_dir = "inputs/"
# Output directory
out_dir = "out/"

mesh = Mesh(in_dir + "mesh.xml")
V_cg = FunctionSpace(mesh, "CG", 1)


### Load model inputs

# Bed
B = Function(V_cg)
File(in_dir + "B.xml") >> B
# Ice thickness
H = Function(V_cg)
File(in_dir + "H.xml") >> H
# Melt
m = Function(V_cg)
File(in_dir + "m.xml") >> m
# Sliding speed
u_b = Function(V_cg)
File(in_dir + "u_b.xml") >> u_b
# Potential at 0 pressure
phi_m = Function(V_cg)
File(in_dir + "phi_m.xml") >> phi_m
# Potential at overburden pressure
phi_0 = Function(V_cg)
File(in_dir + "phi_0.xml") >> phi_0
# Ice overburden pressure
p_i = Function(V_cg)
File(in_dir + "p_i.xml") >> p_i

# Set some parameters for the Newton solver
prm = NonlinearVariationalSolver.default_parameters()
prm['newton_solver']['relaxation_parameter'] = 1.0
prm['newton_solver']['relative_tolerance'] = 1e-6
prm['newton_solver']['absolute_tolerance'] = 1e-6
prm['newton_solver']['error_on_nonconvergence'] = True
prm['newton_solver']['maximum_iterations'] = 50
prm['newton_solver']['linear_solver'] = 'mumps'

# Create the model inputs dictionary
model_inputs = {}
model_inputs['mesh'] = mesh
model_inputs['B'] = B
model_inputs['H'] = H
model_inputs['m'] = m
model_inputs['u_b'] = u_b
model_inputs['phi_m'] = phi_m
model_inputs['phi_0'] = phi_0
model_inputs['p_i'] = p_i
model_inputs['d_bcs'] = []
model_inputs['newton_params'] = prm


### Create the sheet model
sheet_model = SheetModel(model_inputs)

