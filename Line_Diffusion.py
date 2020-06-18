import numpy as np
import numpy.linalg as la
from matplotlib import pyplot as plt
from matplotlib import animation

#Params
v = 0.01
n = 50
dt = 0.01

#min eig of A
lambda_Amin = 0

#builds tridiagonal matrix given center diag and offset diags
def build_tridiag(center, offset, n):
    A = np.zeros((n-1, n-1))
    for i in range(n-1):
        A[i][i] = center[i]
        if (i > 0):
            A[i][i-1] = offset[i]
        if (i < n-2):
            A[i][i+1] = offset[i]
    return A

#ODE EF Solver
def EF(M, init, dt):
    n = len(M) + 1
    ts = [0, 0]
    u = np.zeros((2, n-1))
    u.fill(init)

    for i in range(100):
        u_step = dt*M@u[-1]
        u_next = u[-1]+u_step
        #print(len(u_next), len(u[0]), len(u))
        u = np.append(u, [u_next], axis=0)
        ts.append((i+1)*dt)
    np.delete(u, 0, 0)
    return (ts, u)

#Solves heat Eq
def heat(v, n, dt):
    global lambda_Amin
    #Spacing
    dx = 1/n
    dx2 = dx*dx

    #forms SPD Matrix
    twos = np.zeros(n-1)
    n_ones = np.zeros(n-1)
    twos.fill(2)
    n_ones.fill(-1)
    A = (1/dx2)*build_tridiag(twos, n_ones, n)

    lambda_Amin = min(la.eig(A)[0])
    

    #set up and solve IVP system
    left_mat = -v*A
    return EF(left_mat, 1, dt)

# Animation code from https://jakevdp.github.io/blog/2012/08/18/matplotlib-animation-tutorial/

# First set up the figure, the axis, and the plot element we want to animate
fig = plt.figure()
ax = plt.axes(xlim=(0, 1), ylim=(0, 1.2))
line, = ax.plot([], [], lw=2)

# initialization function: plot the background of each frame
def init():
    line.set_data([], [])
    return line,

# animation function.  This is called sequentially
def animate(i):
    x = np.linspace(1/n, 1-1/n, n-1)
    y = heat(v, n, dt)[1][i]
    line.set_data(x, y)
    return line,

# call the animator.  blit=True means only re-draw the parts that have changed.
anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=100, interval=20, blit=True)
plt.show()

#Setting and testing out over and under max dt
dt_max = (1/n)**2/(2*v)

print("dt_max: " + str(dt_max))

#=======================UNDER dt_max=======================#
dt = 0.95*dt_max

#Print ||u||inf after 100 steps
print("||u100||inf: " + str(max(heat(v, n, dt)[1][-1])))

fig = plt.figure()
ax = plt.axes(xlim=(0, 1), ylim=(0, 1.2))
line, = ax.plot([], [], lw=2)

def init():
    line.set_data([], [])
    return line,

def animate(i):
    x = np.linspace(1/n, 1-1/n, n-1)
    y = heat(v, n, dt)[1][i]
    line.set_data(x, y)
    return line,

anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=100, interval=20, blit=True)
plt.show()

#=======================OVER dt_max=======================#
dt = 1.05*dt_max

#Print ||u||inf after 100 steps
print("||u100||inf: " + str(max(heat(v, n, dt)[1][-1])))

fig = plt.figure()
ax = plt.axes(xlim=(0, 1), ylim=(0, 1.2))
line, = ax.plot([], [], lw=2)

def init():
    line.set_data([], [])
    return line,

def animate(i):
    x = np.linspace(1/n, 1-1/n, n-1)
    y = heat(v, n, dt)[1][i]
    line.set_data(x, y)
    return line,

anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=100, interval=20, blit=True)

plt.show()

#lambdas
lambda_min = np.pi**2
lambda_min_theory = 2*n**2*(1-np.cos(np.pi/n))

#print errors between lambda(A) and above lambdas
print("|lambda_min(A) - lambda~min|: " + str(np.absolute(lambda_Amin - lambda_min)))
print("|lambda_min(A) - lambda_min(theory)|: " + str(np.absolute(lambda_Amin - lambda_min_theory)))
