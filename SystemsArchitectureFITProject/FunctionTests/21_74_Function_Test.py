from SystemsArchitectureFITProject.Functions.Functions import Functions

# Create an instance of the Functions class
functions = Functions()

testLambda_value = 0.955e-8
testTheta_0 = 0.1e-3
testL = 13719818

# Solving for theta_0
result = Functions.function21_74(theta_0=None, lambda_value=testLambda_value, L=testL)
print("Solved θ₀: ", result)

# Solving for lambda_value
result = Functions.function21_74(theta_0=testTheta_0, lambda_value=None, L=testL)
print("Solved λ: ", result)


# Solving for L
result = Functions.function21_74(theta_0=testTheta_0, lambda_value=testLambda_value, L=None)
print("Solved L: ", result)