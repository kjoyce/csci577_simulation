class EulerRichardsonMethod(ODE_integrator)
  def __update_function(f,x,t,dt):
    return (t+dt,x + f(t+dt/2,x+ f(t,x)*dt/2)*dt)
