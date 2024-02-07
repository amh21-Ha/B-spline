# Import numpy for numerical computations
import numpy as np

# Define the system of equations as a function
def f(t, x, y):
    return np.array([y, -np.exp(t-1) - 1])

# Define the initial and final times
t0 = 0
t1 = 1

# Define the number of iterations
n = 30

# Define the step size
h = (t1 - t0) / n

# Define the initial condition for x
x0 = 0

# Define the desired boundary value for x
x1 = 0

# Define the initial interval for s
a = -10
b = 10

# Define the tolerance for the error
tol = 1e-6

# Define a function that solves the initial value problem for a given s
def solve_ivp(s):
    # Initialize the arrays for t, x, and y
    t = np.zeros(n+1)
    x = np.zeros(n+1)
    y = np.zeros(n+1)

    # Set the initial values
    t[0] = t0
    x[0] = x0
    y[0] = s

    # Apply the Runge-Kutta method
    for i in range(n):
        # Compute the intermediate values
        k1 = f(t[i], x[i], y[i])
        k2 = f(t[i] + h/2, x[i] + h*k1[0]/2, y[i] + h*k1[1]/2)
        k3 = f(t[i] + h/2, x[i] + h*k2[0]/2, y[i] + h*k2[1]/2)
        k4 = f(t[i] + h, x[i] + h*k3[0], y[i] + h*k3[1])

        # Update the values of t, x, and y
        t[i+1] = t[i] + h
        x[i+1] = x[i] + h*(k1[0] + 2*k2[0] + 2*k3[0] + k4[0])/6
        y[i+1] = y[i] + h*(k1[1] + 2*k2[1] + 2*k3[1] + k4[1])/6

    # Return the arrays of t, x, and y
    return t, x, y

# Define a function that evaluates the error for a given s
def error(s):
    # Solve the initial value problem for s
    t, x, y = solve_ivp(s)

    # Return the difference between the final value of x and the desired boundary value
    return x[-1] - x1

# Apply the bisection method to find the root of the error function
while abs(b - a) > tol:
    # Compute the midpoint of the interval
    c = (a + b) / 2

    # Evaluate the error function at the endpoints and the midpoint
    fa = error(a)
    fb = error(b)
    fc = error(c)

    # Check if the root is in the left or right subinterval
    if fa * fc < 0:
        # The root is in the left subinterval
        b = c
    elif fc * fb < 0:
        # The root is in the right subinterval
        a = c
    else:
        # The root is either c or an endpoint
        break

# Print the final value of s
print(f"The value of s that satisfies the boundary condition is {c:.6f}")

# Solve the initial value problem for the final value of s
t, x, y = solve_ivp(c)

# Print the solution of the boundary value problem
print(f"The solution of the boundary value problem is x(t) = {x[-1]:.6f} at t = {t[-1]:.6f}")
