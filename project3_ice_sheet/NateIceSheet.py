from dolfin import *
import numpy as np
import pylab as py
from IPython.core.debugger import Tracer
debug_here = Tracer()


###Parameters###

#z_s = 2500                 #surface set at zero
##use for calculateing season variations/firn depth
z_s = 2000
z_b = 0                 #depth of bottom (m)
zboundUp = z_s
zboundLow = z_b
surface_slope = .2E-5   #surface slope of ice
hor_temp_grad = 1E-4    #horizontal temperature gradient
Q_geo = 50E-3           #geothermal heat flow
vert_iceVel_rate = -.2  #vertical ice velocity rate
hor_iceVel_rate = 50.
g = 9.81                #gravity
spy = 31556926          #seconds per year
rho = 911               #ice density
C_p = 2009              #Heat capacity
beta = 9.8E-8           #pressure dependence of melting point
k = 2.1                 #thermal diffusivity of ice
theta_PMP = lambda z: -beta*rho*g*(z_s - z)  #pressure melting point of ice

###Simulation conditionals

save_img = False #save images of steady states
warming = False #apply warming condition after first steady state
killAtSteady = True #stops simulation when steady state is reached
plotLines = False #turn simulation plotting on and off
plotPhase = True #turn phase diagrams on or off

###Time stuff###

##use to see seasonal variations
#timeMax = 30000*spy
#n_steps = 2567
#dt = timeMax/(n_steps)
#t = dt

##use for firn depth
#timeMax = 30000*spy
#n_steps = 30000.
#dt = timeMax/(n_steps*12.)
#t = dt


timeMax = 40000*spy
n_steps = 2000
dt = timeMax/(n_steps)

#timeMax = 30*spy
#n_steps = 30
#dt = timeMax/(n_steps)
#t = dt

#########################Simulation Function###########################

def simulation(zb = z_b, zs = z_s, slope = surface_slope, tempGrad = hor_temp_grad, geo = Q_geo, vertV = vert_iceVel_rate, horV = hor_iceVel_rate):

    z_s = zs
    z_b = zb                #depth of bottom (m)
    zboundUp = z_s
    zboundLow = z_b
    surface_slope = slope   #surface slope of ice
    hor_temp_grad = tempGrad    #horizontal temperature gradient
    Q_geo = geo           #geothermal heat flow
    vert_iceVel_rate = vertV  #vertical ice velocity rate
    hor_iceVel_rate = horV
    
    ###Expressions###

    sigma_expression = '(x[0] - z_b)/(z_s - z_b)'

    sigma = Expression(sigma_expression,
                       z_s = z_s, z_b = z_b)
                       
    u = Expression('u_rate*pow('+sigma_expression+', 4)',
                   z_s = z_s, z_b = z_b, u_rate = hor_iceVel_rate)
                   
    w = Expression(('rate*pow('+sigma_expression+',1)',),
                   z_s = z_s, z_b = z_b, rate = vert_iceVel_rate)
                   
    phi = Expression('-1*rho*g*(z_s - x[0])*200*pow('+sigma_expression+',3)*(1/(z_s - z_b))*slope', rho = rho, g = g, z_s = z_s, z_b = z_b, slope = surface_slope)

    neumannBCCondition = Expression('1-'+sigma_expression,z_s = z_s, z_b = z_b)

    theta_z_s = Expression('base_temp + 5*sin(2*pi*t/spy)', t=0, spy=spy, base_temp = -10)

    ###Mesh###
    mesh = IntervalMesh(1000,z_b,z_s)
    Q = FunctionSpace(mesh, 'CG', 1)

    ###Dirichlet Boundary###

    def dirichletBoundary(x, on_boundary):
        tol = 1E-14
        return on_boundary and (abs(x[0]-z_s) <= tol)
       
    d_BC = DirichletBC(Q, theta_z_s, dirichletBoundary)

    ###Varitional Problem Def###

    T = TrialFunction(Q)
    v = TestFunction(Q)
    theta_last = interpolate(Expression("T0",T0 = -10), Q)
    a = (T*v + k/(rho*C_p)*dt*dot(nabla_grad(T), nabla_grad(v)) + dot(w,nabla_grad(T))*v*dt/spy)*dx
    ##use this for no advection/diffusion
    #a = (T*v + k/(rho*C_p)*dt*dot(nabla_grad(T), nabla_grad(v)))*dx

    ###Matrix assembly###

    A = assemble(a)
    b = None

    T = Function(Q)

    ###Plot setup###
    if plotLines:
        fig = py.figure()
        py.ion()
        ax, = py.plot(theta_last.vector().array(),mesh.coordinates(),'k.')
        py.xlabel("$\\theta$")
        py.ylabel("$z$")
        py.xlim(-10,-7)

    ###Initial Data Setup###
    stateTol = .001
    steadyState = False
    steadyState2 = False
    previousStateVector = py.ones(theta_last.vector().array().shape)*100
    steadyStateTime = "N/A"
    ssTime = -10
    steadyStateTime2 = "N/A"
    update = True
    update2 = True
    temperatures = np.zeros((int(timeMax/dt),len(mesh.coordinates())))
    i=0
    t = dt
    
    ###Return Values###
    depth = mesh.coordinates()[:,0]
    meltPoints = np.zeros((1,n_steps)).flatten()
    
    #########Simulation Loop#########**********************************

    while t <= timeMax:

        ###PDE Solver###
        
        theta_z_s.t = t
        L = (theta_last*v - u*hor_temp_grad*v*dt/spy + phi/(rho*C_p)*v*dt/spy)*dx + (neumannBCCondition*dt/rho/C_p*Q_geo*v)*ds
        ##use this for no advection/diffusion
        #L = (theta_last*v)*dx + (neumannBCCondition*dt/rho/C_p*Q_geo*v)*ds
        b = assemble(L, tensor=b)
        d_BC.apply(A, b)
        solve(A, T.vector(), b)
        t += dt
        theta_last = T
        
        ###Calculations###
        
        npTheta = theta_last.vector().array()
        
        ##check for bad temps
        z = mesh.coordinates()[:,0]
        badIdx = npTheta > theta_PMP(z)
        npTheta[badIdx] = theta_PMP(z)[badIdx]
        theta_last.vector().set_local(npTheta)
        if len(z[badIdx==True]) > 0:
            meltPoints[i] = max(z[badIdx == True])
#        if (t/spy) > 10000:
            #debug_here()
        
        ##steady state check
        
        if(t%(10*dt) == 0):
            currentStateVector = T.vector().array()
            error = py.norm(currentStateVector - previousStateVector)/py.norm(previousStateVector)
            steadyState = error < stateTol
            previousStateVector = currentStateVector
            
        if((t%(10*dt)) == 0 and (update == False) and ((t/spy) > ssTime*1.2)):
            currentStateVector = T.vector().array()
            error = py.norm(currentStateVector - previousStateVector)/py.norm(previousStateVector)
            steadyState2 = error < stateTol
            previousStateVector = currentStateVector
            
        if steadyState and update:
            if warming:
                theta_z_s.base_temp = -8
            ssTime = t/spy
            steadyStateTime = str(t/spy) + " years"
            previousStateVector = py.ones(theta_last.vector().array().shape)*100
            if killAtSteady:
                t = timeMax
            
        if steadyState2 and update2:
            steadyStateTime2 = str(t/spy) + " years"
            
        ###Pull temperature values###
        temperatures[i,:] = npTheta
        i+=1
        
        ###Plotting###
        
        if plotLines:   
            py.title("Time(years): {0}, Max Temp:{1}, Min Temp:{2},\n dt: {3} years, Steady at: {4}, Steady after warming at: {5}".format(t/spy, np.round(max(npTheta), 3), np.round(min(npTheta), 3), (dt/spy), steadyStateTime, steadyStateTime2))
            ax.set_xdata(npTheta)
            py.xlim(-20,0)
            py.ylim(zboundLow,zboundUp)
            py.draw()
    #    if (t/spy)==15:
    #        fig.savefig("test.png", dpi = fig.dpi)
    #        test = False

        ###Images###
        
        if steadyState and update:
            if save_img:
                fig.savefig("FirstSteadyState_{0}_{1}_{2}_{3}_{4}.png".format(surface_slope, hor_temp_grad, Q_geo, vert_iceVel_rate, hor_iceVel_rate))
            update = False
            
        if steadyState2 and update2:
            if save_img:
                fig.savefig("SecondSteadyState_{0}_{1}_{2}_{3}_{4}.png".format(surface_slope, hor_temp_grad, Q_geo, vert_iceVel_rate, hor_iceVel_rate))
            update2 = False
                
    return temperatures, depth, meltPoints, badIdx

if plotPhase:

    ###Horizontal velocity variations
    horVels = np.arange(20,100,5)
    horVels = horVels.tolist()
    meltVectors = list()
    for vel in horVels:
        temps, depths, melts, isMelted = simulation(horV = vel)
        meltVectors.append(isMelted)
        print "tick" 
    print "Ding!"
    fig = py.figure()
    py.contourf(horVels, depths, np.array(meltVectors).T,cmap=py.cm.winter_r)
    py.title("Phase Diagram: Height vs. Horizontal Velocity")
    py.xlabel("Horizontal Ice Flow Velocity")
    py.ylabel("Height from Bed")
    py.ylim(0,50)
    #py.colorbar()
    fig.savefig("horVelPhase.png")
    py.show()
    
    ###Sheet thickness variations
    maxDepth = np.arange(300,3000,300)
    maxDepth = maxDepth.tolist()
    meltVectors = list()
    for dep in maxDepth:
        temps, depths, melts, isMelted = simulation(zs = dep)
        meltVectors.append(isMelted)
        print "tick"
    print "Ding!"
    fig = py.figure()
    py.contourf(maxDepth, depths, np.array(meltVectors).T,cmap=py.cm.winter_r)
    py.title("Phase Diagram: Height vs. Sheet Thickness")
    py.xlabel("Ice Sheet Thickness")
    py.ylabel("Height from Bed")
    py.ylim(0,50)
    #py.colorbar()
    fig.savefig("thicknessPhase.png")
    py.show()
    
    ###Geothermal variations
    geo = np.arange(30E-3,70E-3,5E-3)
    geo = geo.tolist()
    meltVectors = list()
    for G in geo:
        temps, depths, melts, isMelted = simulation(geo = G)
        meltVectors.append(isMelted)
        print "tick"
    print "Ding!"
    fig = py.figure()
    py.contourf(geo, depths, np.array(meltVectors).T,cmap=py.cm.winter_r)
    py.title("Phase Diagram: Height vs. Geothermal Heat Flow")
    py.xlabel("Geothermal Heat Flow")
    py.ylabel("Height from Bed")
    py.ylim(0,50)
    #py.colorbar()
    fig.savefig("QgeoPhase.png")
    py.show()

    ###Vertical velocity variations
    vel = np.arange(-.1,-.5,-0.05)
    vel = vel.tolist()
    meltVectors = list()
    for V in vel:
        temps, depths, melts, isMelted = simulation(vertV = V)
        meltVectors.append(isMelted)
        print "tick"
    print "Ding!"
    fig = py.figure()
    py.contourf(vel, depths, np.array(meltVectors).T,cmap=py.cm.winter_r)
    py.title("Phase Diagram: Height vs. Vertical Ice Velocity")
    py.xlabel("Vertical Ice Velocity")
    py.ylabel("Height from Bed")
    py.ylim(0,50)
    #py.colorbar()
    fig.savefig("vertVelPhase.png")
    py.show()

else:
    simulation()

#temps, depths, melts, isMelted = simulation()
#print "ok"
#print isMelted
#fig = py.figure()
#py.contourf(np.arange(0,timeMax/spy,dt/spy),depths,temps.T,cmap=py.cm.winter)
#py.colorbar()
#py.plot(np.arange(0,timeMax/spy,dt/spy),melts)
#fig.savefig("contourTest.png", )
#py.show()
