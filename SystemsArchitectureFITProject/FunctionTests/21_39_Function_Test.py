from SystemsArchitectureFITProject.Functions.Functions import Functions

# Create an instance of the Functions class
functions = Functions()

# Case 1 Solve for both delta_theta_d_x and delta_theta_d_y
lambda_ = 5
d_x = 6
d_y = 10

print("Solved for both coordinates:", functions.function21_39(lambda_=lambda_, d_x=d_x, d_y=d_y))

# Case 2 Solve for delta_theta_d_x while knowing delta_theta_d_y
lambda_ = 5
d_x = 6
delta_theta_d_y = 10

print("Solved for delta_theta_d_x:", functions.function21_39(delta_theta_d_y=delta_theta_d_y, lambda_=lambda_, d_x=d_x))

# Case 3 Solve for delta_theta_d_y while knowing delta_theta_d_x
lambda_ = 5
d_y = 6
delta_theta_d_x = 10

print("Solved for delta_theta_d_y:", functions.function21_39(delta_theta_d_x=delta_theta_d_x, lambda_=lambda_, d_y=d_y))

# Case 4 Solve for lambda_ while knowing delta_theta_d_x
d_x = 6
delta_theta_d_x = 10

print("Solved for lambda_ with delta_theta_d_x:", functions.function21_39(delta_theta_d_x=delta_theta_d_x, d_x=d_x))

# Case 5 Solve for lambda_ while knowing delta_theta_d_y
d_y = 7
delta_theta_d_y = 10

print("Solved for lambda_ with delta_theta_d_y:", functions.function21_39(delta_theta_d_y=delta_theta_d_y, d_y=d_y))

# Case 6 Solve for d_x
lambda_ = 7
delta_theta_x = 10

print("Solved for d_x:", functions.function21_39(delta_theta_d_x=delta_theta_d_x, lambda_=lambda_))

# Case 7 Solve for d_y
lambda_ = 8
delta_theta_y = 10

print("Solved for d_y:", functions.function21_39(delta_theta_d_y=delta_theta_d_y, lambda_=lambda_))