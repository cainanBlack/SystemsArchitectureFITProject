from SystemsArchitectureFITProject.Functions.Functions import Functions

# Create an instance of the Functions class
functions = Functions()

r_value = functions.calculate_missing_variable(h=1.0, D=10, lamda=3, d_i=2)
print(f"Solved radius (r): {r_value}")

D_value = functions.calculate_missing_variable(h=1.0, r=5, lamda=3, d_i=2)
print(f"Solved diameter (D): {D_value}")

lamda_value = functions.calculate_missing_variable(h=1.0, r=5, D=10, d_i=2)
print(f"Solved wavelength (lambda): {lamda_value}")

d_i_value = functions.calculate_missing_variable(h=1.0, r=5, D=10, lamda=3)
print(f"Solved disk diameter (d_i): {d_i_value}")