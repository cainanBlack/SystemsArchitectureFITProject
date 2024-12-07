from SystemsArchitectureFITProject.Functions.Functions import Functions

# Create an instance of the Functions class
functions = Functions()

# Test Case 1: Solve for r_new
r_new_result = functions.function21_68(r_old=0.2, lambda_new=500e-9, lambda_old=400e-9)
print("Solving for r_new", r_new_result)

# Test Case 2: Solve for r_old
r_old_result = functions.function21_68(r_new=0.25, lambda_new=500e-9, lambda_old=400e-9)
print("Solving for r_old", r_old_result)

# Test Case 3: Solve for lambda_new
lambda_new_result = functions.function21_68(r_new=0.25, r_old=0.2, lambda_old=400e-9)
print("Solving for lambda_new", lambda_new_result)

# Test Case 4: Solve for lambda_old
lambda_old_result = functions.function21_68(r_new=0.25, r_old=0.2, lambda_new=500e-9)
print("Solving for lambda_old", lambda_old_result)


