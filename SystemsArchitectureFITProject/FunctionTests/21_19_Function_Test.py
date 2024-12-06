from SystemsArchitectureFITProject.Functions.Functions import Functions

# Create an instance of the Functions class
functions = Functions()

wp = 1
theta_x = 45

print("Solved for gpf:", functions.function21_19(wp=wp, theta_x=theta_x))

gpf = 1
theta_x = 45

print("Solved for wp:", functions.function21_19(gpf=gpf, theta_x=theta_x))

gpf = 2
wp = 3

print("Solved for theta_x:", functions.function21_19(gpf=gpf, wp=wp))