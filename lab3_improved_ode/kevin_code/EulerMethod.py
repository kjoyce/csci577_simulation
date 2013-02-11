class EulerMethod(ODE_integrator):
  def __update_function(f,x,t,dt):
    return (t+dt,x+f(t,x)*dt)
