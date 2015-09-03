from dolfin import *
from constants import *
from h_solver import *
from phi_solver import *

""" Schoof's simple 1D sheet model with a shocking twist that you won't want
to miss! """

class SheetModel():

  def __init__(self, model_inputs, in_dir = None):

    ### Initialize model inputs

    self.mesh = model_inputs['mesh']
    self.V_cg = FunctionSpace(self.mesh, "CG", 1)
    self.model_inputs = model_inputs
    
    # If an input directory is specified, load model inputs from there. 
    # Otherwise use the specified model inputs dictionary.
    if in_dir:
      model_inputs = self.load_inputs(in_dir)
      
    # Ice thickness    
    self.H = self.model_inputs['H']
    # Bed elevation
    self.B = self.model_inputs['B']
    # Basal sliding speed
    self.u_b = self.model_inputs['u_b']
    # Melt rate
    self.m = self.model_inputs['m']
    # Cavity gap height
    self.h = self.model_inputs['h_init']
    # Potential at 0 pressure
    self.phi_m = self.model_inputs['phi_m']
    # Ice overburden pressure
    self.p_i = self.model_inputs['p_i']
    # Potential at overburden pressure
    self.phi_0 = self.model_inputs['phi_0']
    # Dirichlet boundary conditions
    self.d_bcs = self.model_inputs['d_bcs']
    # Newton solver parameters
    self.newton_params = self.model_inputs['newton_params']
    
    # If there is a dictionary of physical constants specified, use it. 
    # Otherwise use the defaults. 
    if 'constants' in self.model_inputs :
      self.constants = self.model_inputs['constants']
    else :
      self.constants = physical_constants
    

    ### Create some fields

    self.V_cg = FunctionSpace(self.mesh, "CG", 1)

    # Potential
    self.phi = Function(self.V_cg)
    # Effective pressure
    self.N = Function(self.V_cg)
    # Water pressure
    self.p_w = Function(self.V_cg)
    # Pressure as a fraction of overburden
    self.pfo = Function(self.V_cg)
    # Current time
    self.t = 0.0


    ### Create the solver objects

    # Potential solver    
    self.phi_solver = PhiSolver(self)
    # Gap height solver
    self.h_solver = HSolver(self)
    

  # Steps the potential, gap height, and water height forward by dt  
  def solve(self, dt):
    # Step the potential forward by dt with h fixed
    self.phi_solver.solve()
    # Step h forward by dt with phi fixed
    self.h_solver.solve(dt)
    
    
  # Load all model inputs from a directory except for the mesh and initial 
  # conditions on h, h_w, and phi
  def load_inputs(self, in_dir):
    # Bed
    B = Function(self.V_cg)
    File(in_dir + "B.xml") >> B
    # Ice thickness
    H = Function(self.V_cg)
    File(in_dir + "H.xml") >> H
    # Melt
    m = Function(self.V_cg)
    File(in_dir + "m.xml") >> m
    # Sliding speed
    u_b = Function(self.V_cg)
    File(in_dir + "u_b.xml") >> u_b
    # Potential at 0 pressure
    phi_m = Function(self.V_cg)
    File(in_dir + "phi_m.xml") >> phi_m
    # Potential at overburden pressure
    phi_0 = Function(self.V_cg)
    File(in_dir + "phi_0.xml") >> phi_0
    # Ice overburden pressure
    p_i = Function(self.V_cg)
    File(in_dir + "p_i.xml") >> p_i
  
    self.model_inputs['B'] = B
    self.model_inputs['H'] = H
    self.model_inputs['m'] = m
    self.model_inputs['u_b'] = u_b
    self.model_inputs['phi_m'] = phi_m
    self.model_inputs['phi_0'] = phi_0
    self.model_inputs['p_i'] = p_i
    
  
  # Update the effective pressure to reflect current value of phi
  def update_N(self):
    self.N.vector()[:] = self.phi_0.vector().array() - self.phi.vector().array()
    
  
  # Update the water pressure to reflect current value of phi
  def update_pw(self):
    self.p_w.vector()[:] = self.phi.vector().array() - self.phi_m.vector().array()
    
  
  # Update the pfo to reflect the current value of phi
  def update_pfo(self):
    # Update water pressure
    self.update_pw()
    # Compute overburden pressure
    self.pfo.vector()[:] = self.p_w.vector().array() / self.p_i.vector().array()
    
  
  # Updates all fields derived from phi
  def update_phi(self):
    self.update_N()
    self.update_pfo()
    