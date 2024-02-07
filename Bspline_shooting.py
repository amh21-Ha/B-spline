#!/usr/bin/python3
import numpy as np

# Define the system of equations and the boundary conditions
def f(t, x, y):
    return y

def g(t, x, y):
    return x + np.exp(t-1) + 1

def x0(s):
    return 0

def y0(s):
    return s

def x1(s):
    return 0

# Define the interval and the tolerance for the bisection method
a = -10
b = 10
tol = 1e-6

# Define the step size and the number of steps for the Runge-Kutta method
h = 0.1
n = int(1/h)

# Define a function that solves the IVP for a given value of s and returns x(1)
def solve_ivp(s):
    # Initialize the arrays for t, x, and y
    t = np.zeros(n+1)
    x = np.zeros(n+1)
    y = np.zeros(n+1)
    # Set the initial conditions
    t[0] = 0
    x[0] = x0(s)
    y[0] = y0(s)
    # Apply the Runge-Kutta method
    for i in range(n):
        # Compute the intermediate values
        k1 = h * f(t[i], x[i], y[i])
        l1 = h * g(t[i], x[i], y[i])
        k2 = h * f(t[i] + h/2, x[i] + k1/2, y[i] + l1/2)
        l2 = h * g(t[i] + h/2, x[i] + k1/2, y[i] + l1/2)
        k3 = h * f(t[i] + h/2, x[i] + k2/2, y[i] + l2/2)
        l3 = h * g(t[i] + h/2, x[i] + k2/2, y[i] + l2/2)
        k4 = h * f(t[i] + h, x[i] + k3, y[i] + l3)
        l4 = h * g(t[i] + h, x[i] + k3, y[i] + l3)
        # Update the values of t, x, and y
        t[i+1] = t[i] + h
        x[i+1] = x[i] + (k1 + 2*k2 + 2*k3 + k4)/6
        y[i+1] = y[i] + (l1 + 2*l2 + 2*l3 + l4)/6
    # Return the value of x(1)
    return x[n]

# Apply the bisection method to find the root of x(1) - 0
while abs(b - a) > tol:
    # Compute the midpoint and the value of x(1) at the midpoint
    c = (a + b)/2
    xc = solve_ivp(c)
    # Check the sign of x(1) and update the interval accordingly
    if xc > 0:
        b = c
    else:
        a = c

# Print the final value of s and x(1)
print("The value of s is", c)
print("The value of x(1) is", xc)
