from dolfin import *
from constants import *
from gap_solver import *
from potential_solver import *

""" Schoof's simple 1D sheet model with a shocking twist that you won't want
to miss! """

class SheetModel():

  def __init__(self, model_inputs):

    ### Get some fields

    # Mesh
    self.mesh = model_inputs['mesh']
    # Ice thickness    
    self.H = model_inputs['H']
    # Bed elevation
    self.B = model_inputs['B']
    # Basal sliding speed
    self.u_b = model_inputs['u_b']
    # Melt rate
    self.m = model_inputs['m']
    # Potential at 0 pressure
    self.phi_m = model_inputs['phi_m']
    # Ice overburden pressure
    self.p_i = model_inputs['p_i']
    # Potential at overburden pressure
    self.phi_0 = model_inputs['phi_0']
    # Dirichlet boundary conditions
    self.d_bcs = model_inputs['d_bcs']
    # Newton solver parameters
    self.newton_params = model_inputs['newton_params']
    
    
    ### Create some fields
    
    self.V_cg = FunctionSpace(self.mesh, "CG", 1)
    # Cavity gap height
    self.h = Function(self.V_cg)
    # Sheet water height
    self.h_w = Function(self.V_cg)
    # Potential
    self.phi = Function(self.V_cg)
    # Effective pressure
    self.N = Function(self.V_cg)
    # Current time
    self.t = 0.0


    ### Create the solver objects

    # Potential solver    
    self.phi_solver = PotentialSolver(self)
    # Gap height solver
    self.gap_solver = GapSolver(self)
    

  # Steps the potential, gap height, and water height forward by dt  
  def solve(self, dt):
    # Step h forward by dt with phi fixed
    self.gap_solver.solve(dt)
    # Step the potential forward by dt with h fixed
    self.phi_solver.solve(dt)