# B-spline-m B-spline numerical solution with python
The B-spline method is a numerical technique that uses cubic B-spline basis functions to approximate the solution of a boundary value problem. The idea is to divide the interval [0, 1] into n subintervals and construct a piecewise cubic polynomial that satisfies the boundary conditions and the differential equation at each node. The coefficients of the polynomial are determined by solving a system of linear equations.

To solve the problem, I will use the following steps:

Step 1: Define the B-spline basis functions and their derivatives on each subinterval.

Step 2: Construct the collocation matrix and the right-hand side vector using the B-spline basis functions and the differential equation.

Step 3: Solve the linear system for the coefficients of the B-spline approximation.

Step 4: Evaluate the B-spline approximation at any desired point in the interval [0, 1].
