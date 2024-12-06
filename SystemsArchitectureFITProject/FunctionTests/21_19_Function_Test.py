from SystemsArchitectureFITProject.Functions.Functions import Functions

# Create an instance of the Functions class
functions = Functions()

# Case 1 Solve for GPF
wp = 1
theta_x = 45

print("Solved for gpf:", functions.function21_19(wp=wp, theta_x=theta_x))

# Case 2 Solve for WP
gpf = 1
theta_x = 45

print("Solved for wp:", functions.function21_19(gpf=gpf, theta_x=theta_x))

# Case 3 Solve for Theta_X
gpf = 2
wp = 3

print("Solved for theta_x:", functions.function21_19(gpf=gpf, wp=wp))