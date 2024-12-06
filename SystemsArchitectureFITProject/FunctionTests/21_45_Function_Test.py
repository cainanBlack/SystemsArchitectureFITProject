from SystemsArchitectureFITProject.Functions.Functions import Functions

# Create an instance of the Functions class
functions = Functions()

# Case 1: Solve for Φ_n^K
Phi_result = functions.function21_45(missing='Phi', C_n=1e-14, K=1.0)

# Case 2: Solve for C_n
C_n_result = functions.function21_45(missing='C_n', Phi=3.3e-16, K=1.0)

# Case 3: Solve for K
K_result = functions.function21_45(missing='K', Phi=3.3e-16, C_n=1e-14)

# Case 4: Calculate when no variable is missing
direct_result = functions.function21_45(missing=None, Phi=3.3e-16, C_n=1e-14, K=1.0)
