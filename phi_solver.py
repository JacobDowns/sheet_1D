from dolfin import *
from constants import *

""" The potential solver solves for the potential phi. """

class PhiSolver():

  def __init__(self, model):
    
    ### Get a few fields from the model    
    
    # Get melt rate
    m = model.m
    # Sheet height
    h = model.h
    # Basal sliding speed
    u_b = model.u_b
    # Potential
    phi = model.phi
    # Potential at overburden pressure
    phi_0 = model.phi_0
       
    
    ### Set up the model
  
    V_cg = model.V_cg     
    # Regularization parameter
    phi_reg = Constant(1e-16)
    # Effective pressure expression
    N = phi_0 - phi
    # Potential gradient
    dphi_ds = phi.dx(0)
    # Flux
    q = -Constant(k) * h**alpha * abs(dphi_ds + phi_reg)**delta * dphi_ds
    # Opening term
    w = u_b * (Constant(h_r) - h) / Constant(l_r)
    # Closing term
    v = Constant(A) * h * N**3
    
    
    ### Set up the variational problem    
    
    # Trial function
    dphi = TrialFunction(V_cg)
    # Test function
    theta = TestFunction(V_cg)
    # Variational form 
    F = (-theta.dx(0) * q + (w - v - m) * theta) * dx
    # Jacobian
    J = derivative(F, phi, dphi)
    
    
    ### Assign local variables
    
    self.F = F
    self.J = J
    self.model = model


  # Solve for the potential
  def solve(self):
    # Solve for the potential 
    solve(self.F == 0, self.model.phi, self.model.d_bcs, J = self.J, solver_parameters = self.model.newton_params)
    # Update fields derived from phi
    self.model.update_phi()
