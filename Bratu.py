import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt

#spacing and constant params
n = 15
sigma = 2.5

#logs finite diff approximations of u'(0)
log = []

#logs sequential differences between approximations of doubling n
errors = []

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

#solves bratu eq
def bratu(n, sigma):
    #spacing
    h = 1/n
    h2i = 1/(h*h)

    #offset diags
    b = np.zeros(n-1)
    for i in range(n-1):
        b[i] = 1
    a = -2*b
    
    #setup A
    A = h2i*build_tridiag(a, b, n)

    #setup J
    x = h*np.array(range(1, n))
    c = -2*b + sigma*h*h*np.exp(x)
    J = h2i*build_tridiag(c, b, n)

    #setup u
    u = b*0
    
    #Newton Iter
    nf = 1
    k = 0
    while (nf > 10**-9):
        f = A@u + sigma*np.exp(u)
        c = -2*b + sigma*h*h*np.exp(u)
        J = h2i*build_tridiag(c, b, n)
        s = la.solve(J, f)
        u = u-s
        nf = la.norm(f, 2)
        k += 1

    #plot curve, update log, and record iters
    plt.plot(x, u)
    plt.title("Bratu Solution, u vs x")
    log.append(u[0]/h)
    print("It takes %s iterations to get a residual norm less than 10^9" %(k))
    return u

#solve curve with n and sigma
bratu(n, sigma)
plt.show()

#solve sequentially doubled nodes
bratu(2*n, sigma)
bratu(4*n, sigma)
bratu(8*n, sigma)
bratu(16*n, sigma)
bratu(32*n, sigma)

plt.show()

#Print approximations for u'(0)
print('\nFor n=15 and sigma=2.5,')
m = 1
for i in log:
    print(str(m)+"n u'(0): " + str(i))
    m = m*2
print('')

#Collect errors
for i in range(1, len(log)):
    errors.append(np.absolute(log[i]-log[i-1]))

#Plot errors, gets straight line
hs = [1/15, 1/30, 1/60, 1/120, 1/240, 1/480]
plt.plot(hs[1:], errors)
plt.title("Convergence of u'(0) vs Spacing")
plt.show()

#Richardson extrapolation for u'(0) approximations
def richardson(index):
    return 2*log[i]-log[i-1]

print('\nNow with Richardson Extrapolation,')

#collect new errors with richardson extrapolation
richardson_errs = []
prev = 0
curr = 0
for i in range(1, len(log)):
    curr = richardson(i)
    print(str(2**(i))+"n: " + str(richardson(i)))
    if (i > 1):
        richardson_errs.append(np.absolute(curr-prev))
    prev = richardson(i)

#plot these new errors, is parabolic
plt.plot(hs[2:], richardson_errs)
plt.title("Convergence of u'(0) with Richardson Extrapolation vs Spacing")
plt.show()

sigma = 3.507
n = 17
bratu(n, sigma)
plt.show()
