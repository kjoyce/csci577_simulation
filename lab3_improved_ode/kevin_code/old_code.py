for dt in dt_range:
  rk =  oi.RungeKutta     (f,y0,a,b,dt)
  em =  oi.Euler          (f,y0,a,b,dt)
  erm = oi.EulerRichardson(f,y0,a,b,dt)
  figure()
  errors = []
  plot_sym = ['+','x','.']
  for method in (em,erm,rk):
    (t,y) = method.integrate()
    plot(t,y[0],plot_sym.pop())
    energy = oscillator_energy(y[0],y[1],k,m)
    errors.append(abs(energy[-1] - energy[0])/energy[0])
  plot(arange(a,b,dt),oscillator_analytic(arange(a,b,dt),k,y0[0],m).T[0],'k-')
  legend((r'Euler, Error={:.4f}'.format(errors[0]),
          r'Euler-Richardson, Error = {:.4f}'.format(errors[1]),
	  r'Runge-Kutta, Error = {:.4f}'.format(errors[2]),'Analytic'),loc='lower left')
  show()




