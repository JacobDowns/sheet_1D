from dolfin import *
from constants import *

""" The potential solver solves for the potential phi. """

class PotentialSolver():

  def __init__(self, model):
    
    ### Get a few fields from the model    
    
    # Get melt rate
    m = model.m
    # Sheet water height
    h_w = model.h_w
    # Potential
    phi = model.phi
  
    ### Create a few additional fields
  
    V_cg = model.V_cg     
    # Sheet water at previous time step
    h_w0 = Function(V_cg)
    # Time step
    dt = Constant(1.0)
    # Potential gradient
    dphi_ds = phi.dx(0)
    # Flux
    q = k * h_w**alpha * abs(dphi_ds)**delta * dphi_ds
    
    
    ### Set up the variational problem    
    
    # Trial function
    dh_w = TrialFunction(V_cg)
    # Test function
    w = TestFunction(V_cg)
    # Variational form 
    F = ((h_w0 - h_w) * w + dt * q * w.dx(0) - dt * m * w) * dx
    # Jacobian
    J = derivative(F, h_w, dh_w)
    
    
    ### Assign local variables
    
    self.dt = dt
    self.q = q
    self.F = F
    self.J = J

  # Solve for the potential
  def solve(self, dt):
    
    # Assign the time step
    self.dt.assign(dt)
    # Solve for the potential 
    solve(self.F == 0, model.phi, model.d_bcs, J = self.J, solver_parameters = model.newton_params)
    # Update the effective pressure 
    model.N.vector()[:] = model.phi0.vector().array() - model.phi.vector().array()