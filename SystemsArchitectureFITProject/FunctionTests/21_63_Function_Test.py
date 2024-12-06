from SystemsArchitectureFITProject.Functions.Functions import Functions

functions = Functions()

# Example 1: Solve for Gamma_p given Delta_x and r0
Delta_x = 2.0
r0 = 1.0
Gamma_p_solution = functions.function21_63(Delta_x=Delta_x, r0=r0)
print(f"Solution for Gamma_p: {Gamma_p_solution}")

# Example 2: Solve for Delta_x given Gamma_p and r0
Gamma_p = 0.5
Delta_x_solution = functions.function21_63(Gamma_p=Gamma_p, r0=1.0)
print(f"Solution for Delta_x: {Delta_x_solution}")

# Example 3: Solve for r0 given Gamma_p and Delta_x
r0_solution = functions.function21_63(Gamma_p=Gamma_p, Delta_x=2.0)
print(f"Solution for r0: {r0_solution}")