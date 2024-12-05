from SystemsArchitectureFITProject.Functions.Functions import Functions

# Create an instance of the Functions class
functions = Functions()

# Case 1 Solve for both missing coordinates with every other value given
lambda_ = 5
z = 6
d = 10

print("Solved for both coordinates:", functions.function21_40(lambda_=lambda_, z=z, d=d))

# Case 2 Solve for delta_r_d while knowing delta_theta_d
lambda_ = 5
z = 6
d = 10
delta_theta_d = 3

print("Solved for delta_r_d:", functions.function21_40(delta_theta_d=delta_theta_d, lambda_=lambda_, z=z, d=d))

# Case 3 Solve for delta_theta_d while knowing delta_r_d
lambda_ = 5
d = 10
delta_r_d = 3

print("Solved for delta_theta_d:", functions.function21_40(delta_r_d=delta_r_d, lambda_=lambda_, d=d))

# Case 4 Solve for d while knowing delta_r_d
lambda_ = 5
z = 6
delta_r_d = 10

print("Solved for d with delta_r_d:", functions.function21_40(delta_r_d=delta_r_d, lambda_=lambda_, z=z))

# Case 5 Solve for d while knowing delta_theta_d
lambda_ = 5
delta_theta_d = 10

print("Solved for d with delta_theta_d:", functions.function21_40(delta_theta_d=delta_theta_d, lambda_=lambda_))

# Case 6 Solve for lambda while knowing delta_r_d
z = 6
d = 10
delta_r_d = 3

print("Solved for lambda with delta_r_d:", functions.function21_40(delta_r_d=delta_r_d, z=z, d=d))

# Case 7 Solve for lambda while knowing delta_theta_d
d = 10
delta_theta_d = 3

print("Solved for lambda with delta_theta_d:", functions.function21_40(delta_theta_d=delta_theta_d, d=d))

# Case 8 Solve for z
lambda_ = 5
delta_r_d = 3
d = 10
delta_theta_d = 6

print("Solved for z:", functions.function21_40(delta_r_d=delta_r_d, delta_theta_d=delta_theta_d, lambda_=lambda_, d=d))