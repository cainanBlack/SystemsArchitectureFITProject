from SystemsArchitectureFITProject.Functions.Functions import Functions

# Create an instance of the Functions class
functions = Functions()

# Test Case 1: Solve for Phi
Phi_result = functions.function21_4(missing='Phi', C_n=1e-14, K=1.0, L_0=100.0, K_m=10.0)

# Test Case 2: Solve for C_n
C_n_result = functions.function21_4(missing='C_n', Phi=3.24e-30, K=1.0, L_0=100.0, K_m=10.0)

# Test Case 3: Solve for L_0
L_0_result = functions.function21_4(missing='L_0', Phi=3.24e-30, C_n=1e-14, K=1.0, K_m=10.0)

# Test Case 4: Solve for K_m
K_m_result = functions.function21_4(missing='K_m', Phi=3.24e-30, C_n=1e-14, K=1.0, L_0=100.0)

# Results
print(f"Φ_n^V(K) (solving for Phi): {Phi_result:.4e}")
print(f"C_n (solving for C_n): {C_n_result:.4e}")
print(f"L_0 (solving for L_0): {L_0_result:.4e}")
print(f"K_m (solving for K_m): {K_m_result:.4e}")
