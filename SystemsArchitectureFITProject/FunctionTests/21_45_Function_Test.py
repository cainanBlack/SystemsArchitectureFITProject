from SystemsArchitectureFITProject.Functions.Functions import Functions

# Create an instance of the Functions class
functions = Functions()

# Case 1: Solve for Φ_n^K
Phi_result = functions.function21_45(C_n=1e-14, K=1.0)
print("Solved for missing Phi:", Phi_result)

# Case 2: Solve for C_n
C_n_result = functions.function21_45(Phi=3.3e-16, K=1.0)
print("Solved for missing C_n:", C_n_result)

# Case 3: Solve for K
K_result = functions.function21_45(Phi=3.3e-16, C_n=1e-14)
print("Solved for missing K:", K_result)


