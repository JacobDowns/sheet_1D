"""
1D Sheet model sim
"""

from dolfin import *
from sheet_model import *
from constants import *
from pylab import *

# Model input directory
in_dir = "inputs_slope/"
# Output directory
out_dir = "out/"

mesh = Mesh(in_dir + "mesh.xml")
V_cg = FunctionSpace(mesh, "CG", 1)

### initial conditions

# Initial sheet height
h_init = Function(V_cg)
h_init.interpolate(Constant(0.05))

# Set some parameters for the Newton solver
prm = NonlinearVariationalSolver.default_parameters()
prm['newton_solver']['relaxation_parameter'] = 1.0
prm['newton_solver']['relative_tolerance'] = 1e-6
prm['newton_solver']['absolute_tolerance'] = 1e-6
prm['newton_solver']['error_on_nonconvergence'] = False
prm['newton_solver']['maximum_iterations'] = 50
prm['newton_solver']['linear_solver'] = 'mumps'

# 0 pressure bc at margin
def margin(x, on_boundary):
  return on_boundary and near(x[0], 0.0)
  
bc = DirichletBC(V_cg, 0.0, margin)

model_inputs = {}
model_inputs['mesh'] = mesh
model_inputs['h_init'] = h_init
model_inputs['d_bcs'] = [bc]
model_inputs['newton_params'] = prm
model_inputs['constants'] = physical_constants

spd = physical_constants['spd']
dt = 60.0 * 30.0
T = 25.0 * spd

sheet_model = SheetModel(model_inputs, in_dir)

xs = project(Expression("x[0]"), V_cg)
indexes = xs.vector().array().argsort()
xs = xs.vector().array()[indexes] / 1000.0



while sheet_model.t < T:
  sheet_model.solve(dt)
  
  #plot(sheet_model.pfo, interactive = False)
  #plot(sheet_model.h, interactive = False)
  
hs = sheet_model.h.vector().array()[indexes]
pfos = sheet_model.pfo.vector().array()[indexes]

plot(xs, phi_ms, 'b')
plot(xs, phi_0s, 'r')
plot(xs, phis, 'k')

show()

z = zeros(len(xs))
title("Pressure as a Fraction of Overburden (k = 0.005)")
plot(xs, pfos, 'k')
plot(xs, z, 'r')
xlabel("x (km)")
ylabel("pfo")
show()
