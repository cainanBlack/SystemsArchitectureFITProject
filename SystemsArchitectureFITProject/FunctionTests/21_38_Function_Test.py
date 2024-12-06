from SystemsArchitectureFITProject.Functions.Functions import Functions

# Create an instance of the Functions class
functions = Functions()

# Case 1 Solve for both missing coordinates with every other value given
lambda_ = 5
z = 10
d_x = 7
d_y = 8

print("Solved for both coordinates:", functions.function21_38(lambda_=lambda_, z=z, d_x=d_x, d_y=d_y))

# Case 2 Solve for delta_x while knowing delta_y
lambda_ = 5
z = 12
d_x = 6
delta_y = 6.25

print("Solved for delta_x:", functions.function21_38(delta_y=delta_y, lambda_=lambda_, z=z, d_x=d_x, d_y=d_y))

# Case 3 Solve for delta_y while knowing delta_x
lambda_ = 5
z = 12
d_y = 6
delta_x = 6.25

print("Solved for delta_y:", functions.function21_38(delta_x=delta_x, lambda_=lambda_, z=z, d_x=d_x, d_y=d_y))

# Case 4 Solve for lambda_ while knowing delta_x
z = 12
d_x = 7
delta_x = 6.25

print("Solved for lambda from delta_x:", functions.function21_38(delta_x=delta_x, z=z, d_x=d_x))

# Case 5 Solve for lambda_ while knowing delta_y
z = 13
d_y = 7
delta_y = 6.25

print("Solved for lambda from delta_y:", functions.function21_38(delta_y=delta_y, z=z, d_y=d_y))

# Case 6 Solve for z while knowing delta_x
delta_x=5
lambda_=8
d_x=10

print("Solved for z from delta_x:", functions.function21_38(delta_x=delta_x, lambda_=lambda_, d_x=d_x))

# Case 7 Solve for z while knowing delta_y
delta_y=5
lambda_=7
d_y=10

print("Solved for z from delta_y:", functions.function21_38(delta_y=delta_y, lambda_=lambda_, d_y=d_y))

# Case 8 Solve for d_x
lambda_ = 5
z = 8
delta_x = 10

print("Solved for d_x:", functions.function21_38(delta_x=delta_x, lambda_=lambda_, z=z))

# Case 9 Solve for d_y
lambda_ = 6
z = 8
delta_y = 10

print("Solved for d_y:", functions.function21_38(delta_y=delta_y, lambda_=lambda_, z=z))