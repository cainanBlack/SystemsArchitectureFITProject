from SystemsArchitectureFITProject.Functions.Functions import Functions

# Create an instance of the Functions class
functions = Functions()

# Test Case 1: Solving for missing phi
rho = 1.0
theta = 1.0
a = [2, 3, None]  # Missing the third coefficient (a[2] is None)
Z = [lambda rho, theta: 2*rho, lambda rho, theta: 3*theta, lambda rho, theta: rho*theta]

# Solve for missing phi
missingPhi = functions.function21_85(rho=rho, theta=theta, a=a, Z=Z)
print("Solved for missing phi:", missingPhi)

# Test Case 2: Solving for missing rho
phi = 10
theta = 1.0
a = [2, 3, 4]
Z = [lambda rho, theta: 2*rho, lambda rho, theta: 3*theta, lambda rho, theta: rho*theta]

missingRho = functions.function21_85(phi=phi, theta=theta, a=a, Z=Z)
print("Solved for missing rho:", missingRho)

# Test Case 3: Solving for missing theta
phi = 10
rho = 1.0
a = [2, 3, 4]
Z = [lambda rho, theta: 2*rho, lambda rho, theta: 3*theta, lambda rho, theta: rho*theta]

missingTheta = functions.function21_85(phi=phi, rho=rho, a=a, Z=Z)
print("Solved for missing theta:", missingTheta)

# Test Case 4: Solving for missing a_i coefficients
phi = 10
rho = 1.0
theta = 1.0
a = [2, None, 4]  # Missing the second coefficient (a[1] is None)
Z = [lambda rho, theta: 2*rho, lambda rho, theta: 3*theta, lambda rho, theta: rho*theta]

missingA = functions.function21_85(phi=phi, rho=rho, theta=theta, Z=Z, a=a)
print("Solved for missing a:", missingA)