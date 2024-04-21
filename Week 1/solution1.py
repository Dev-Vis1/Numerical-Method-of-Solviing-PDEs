"""
from numpy import arange
def ftcsadvect(x,dt,T,rho0):
#-- compute spatial step size
    dx= x[1]-x[0]

    #-- initial condition; set up future time array
    u0=rho0.copy()
    u1=u0.copy()
    nu= dt/(2.*dx)

    #-- now the time loop
    for t in arange(0,T+dt,dt):
        u1[1:(len(x)-1)]= u0[1:(len(x)-1)] - nu* (u0[2:len(x)]-u0[0:(len(x)-2)])
        u1[0] = u0[0] - nu* (u0[1] - u0[len(x)-1])
        u1[len(x)-1]= u0[(len(x)-1)] - nu* (u0[0]-u0[len(x)-2])
        u0=u1.copy()

    #-- set output value
    rho = u1  
    return rho

# initial condition
from numpy import pi, exp, cos
def init1(x):
    u0=exp(-(20*x)**2/2)
    return u0
def init2(x):
    sigm= 0.1; Ka= pi/sigm
    u0= exp(-(x**2)/(2*sigm**2)) * cos(Ka*x)
    return u0

def testlf():
    from numpy import arange, pi, exp, cos, inf
    from numpy.linalg import norm
    import matplotlib.pyplot as plt
    
    #-- initial values/settings
    dx= 1./100.              #- spatial step size
    courant= 0.5             #- Courant No.
    dt= courant* dx          #- time step size
    Tend= .5                 #- final time
    x= arange(-1,1+dx,dx)    #- spatial grid
    
    #-- initial conditions
    u0 = init1(x)
    
    #-- Lax-Friedrichs advection scheme
    uftcs = ftcsadvect(x,dt,Tend,u0)
    
    #-- compute exact solution to compare with
    uexact = init1(x-Tend)
    
    #-- plot the result, the initial condition, and the exact solution
    fig = plt.figure(1)
    h1=plt.plot(x,u0,linewidth=2, c=[0.7, 0.7, 0.7], label='init. cond.')
    h2=plt.plot(x,uexact,linewidth=2, c=[0.5, 0.7, 0.4], label='exact sol.')
    h3=plt.plot(x,uftcs,linewidth=2, c='red', label='numeric sol.')
    plt.legend(loc='upper left')
    plt.title('Linear Advection, DX=' + '%6.4f' % (dx) + ', DT=' + '%6.4f' % (dt))
    plt.xlabel('x')
    plt.ylabel('u')
  
    #-- compute error norms:
    infftcs= norm((uexact-uftcs),inf)/norm(uexact,inf)
    print('*** INFO: relative error in inf-norm ***')
    print('          FTCS method: ' + '%6.4f' % (infftcs))
    twoftcs= norm((uexact-uftcs))/norm(uexact)
    print('*** INFO: relative error in two-norm ***')
    print('          FTCS method: ' + '%6.4f' % (twoftcs))
    
    plt.show()
    if __name__=="__main__":
    testlf()
"""


from numpy import arange, exp, inf
import matplotlib.pyplot as plt
from numpy.linalg import norm

def leapfrogadvect(x, dt, T, rho0, a=0.05):
    dx = x[1] - x[0]
    N = len(x)
    u0 = rho0.copy()
    nu = dt / (2.0 * dx)
    u1 = u0.copy()  # Initialize rho at half step
    
    for t in arange(0, T + dt, dt):
        u1[1:(N - 1)] = u0[1:(N - 1)] - nu * (u0[2:N] - u0[0:(N - 2)])
        u1[0] = u0[0] - nu * (u0[1] - u0[N - 1])
        u1[N - 1] = u0[N - 1] - nu * (u0[0] - u0[N - 2])
        
        # Apply Asselin Filter
        u0[1:(N - 1)] = u1[1:(N - 1)] + a * (u0[2:N] - 2 * u0[1:(N - 1)] + u0[0:(N - 2)])
        u0[0] = u1[0] + a * (u0[1] - 2 * u0[0] + u0[N - 1])
        u0[N - 1] = u1[N - 1] + a * (u0[0] - 2 * u0[N - 1] + u0[N - 2])
    
    return u0

# initial condition
def init1(x):
    u0 = exp(-(20 * x) ** 2 / 2)
    return u0

def test_leapfrog():
    dx = 1. / 100.  # spatial step size
    courant = 0.5  # Courant number
    dt = courant * dx  # time step size
    Tend = 0.5  # final time
    x = arange(-1, 1 + dx, dx)  # spatial grid
    
    # initial conditions
    u0 = init1(x)
    
    # Leap-Frog advection scheme with Asselin Filter
    rho = leapfrogadvect(x, dt, Tend, u0, a=0.05)
    
    # compute exact solution to compare with
    uexact = init1(x - Tend)
    
    # plot the result, the initial condition, and the exact solution
    fig = plt.figure(1)
    h1 = plt.plot(x, u0, linewidth=2, c=[0.7, 0.7, 0.7], label='init. cond.')
    h2 = plt.plot(x, uexact, linewidth=2, c=[0.5, 0.7, 0.4], label='exact sol.')
    h3 = plt.plot(x, rho, linewidth=2, c='red', label='numeric sol.')
    plt.legend(loc='upper left')
    plt.title('Linear Advection, DX=' + '%6.4f' % dx + ', DT=' + '%6.4f' % dt)
    plt.xlabel('x')
    plt.ylabel('u')
    
    # compute error norms
    inff = norm((uexact - rho), inf) / norm(uexact, inf)
    print('*** INFO: relative error in inf-norm ***')
    print('          Leap-Frog method with Asselin Filter: ' + '%6.4f' % inff)
    twof = norm((uexact - rho)) / norm(uexact)
    print('*** INFO: relative error in two-norm ***')
    print('          Leap-Frog method with Asselin Filter: ' + '%6.4f' % twof)
    
    plt.show()

if __name__ == "__main__":
    test_leapfrog()