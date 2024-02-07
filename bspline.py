#!/usr/bin/python3
# Import the libraries
import numpy as np
import matplotlib.pyplot as plt

# Define the parameters
n = 10 # number of subintervals
h = 1/n # length of subinterval
t = np.linspace(0, 1, n+1) # nodes
f = lambda t: -np.exp(t-1) - 1 # right-hand side function

# Define the B-spline basis functions and their derivatives
def N(i, t):
    # B-spline basis function of degree 3
    if i == 0:
        return (1-t)**3/6
    elif i == 1:
        return (3*t**3 - 6*t**2 + 4)/6
    elif i == 2:
        return (-3*t**3 + 3*t**2 + 3*t + 1)/6
    elif i == 3:
        return t**3/6
    else:
        return 0

def dN(i, t):
    # Derivative of B-spline basis function of degree 3
    if i == 0:
        return -0.5*(1-t)**2
    elif i == 1:
        return 1.5*t**2 - 2*t
    elif i == 2:
        return -1.5*t**2 + t + 0.5
    elif i == 3:
        return 0.5*t**2
    else:
        return 0

def ddN(i, t):
    # Second derivative of B-spline basis function of degree 3
    if i == 0:
        return 1 - t
    elif i == 1:
        return 3*t - 2
    elif i == 2:
        return -3*t + 1
    elif i == 3:
        return t
    else:
        return 0

# Construct the collocation matrix and the right-hand side vector
A = np.zeros((n+1, n+1)) # collocation matrix
b = np.zeros(n+1) # right-hand side vector
for i in range(n+1):
    for j in range(n+1):
        A[i, j] = ddN(j, t[i]) - dN(j, t[i]) # coefficient of x_j at node t_i
    b[i] = f(t[i]) # right-hand side at node t_i
# Apply the boundary conditions
A[0, 0] = 1 # x_0 = 0
A[0, 1:] = 0 # x_0 = 0
b[0] = 0 # x_0 = 0
A[-1, -1] = 1 # x_n = 0
A[-1, :-1] = 0 # x_n = 0
b[-1] = 0 # x_n = 0

# Solve the linear system for the coefficients
x = np.linalg.solve(A, b)

# Evaluate the B-spline approximation at any desired point
def S(t):
    # B-spline approximation of degree 3
    s = 0
    for i in range(n+1):
        s += x[i] * N(i, t)
    return s

# Plot the B-spline approximation and the exact solution
t_plot = np.linspace(0, 1, 100) # points for plotting
S_plot = np.vectorize(S)(t_plot) # B-spline approximation at plot points
exact = lambda t: -np.exp(t-1) # exact solution
exact_plot = exact(t_plot) # exact solution at plot points
plt.plot(t_plot, S_plot, label='B-spline approximation')
plt.plot(t_plot, exact_plot, label='Exact solution')
plt.scatter(t, x, label='B-spline coefficients')
plt.xlabel('t')
plt.ylabel('x')
plt.legend()
plt.show()
