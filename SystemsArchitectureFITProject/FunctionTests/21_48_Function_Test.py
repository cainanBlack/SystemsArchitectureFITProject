from SystemsArchitectureFITProject.Functions.Functions import Functions

# Create an instance of the Functions class
functions = Functions()

# Test Case 1: Solve for Phi
print("Solve for Phi:", functions.function21_48(C_n=1e-14, K=1.0, K_0=100.0, K_m=10.0))

# Test Case 2: Solve for C_n
C_n_result = functions.function21_48(Phi=3.24e-30, K=1.0, K_0=100.0, K_m=10.0)
print("Solve for C_n:", C_n_result)