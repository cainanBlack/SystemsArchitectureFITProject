from SystemsArchitectureFITProject.Functions.Functions import Functions

# Create an instance of the Functions class
functions = Functions()

# Values for Test
testPhi = 1               # Simple example value for Φ_n(k)
testGamma = 5 + 2j        # Γₙ(r)
testK = (1.0, 1.0, 1.0)   # Example value for k (wave vector in 3D space)

# Case 1: Given values and solve for the missing gamma_n
missing_gamma_n = functions.function21_59(phi=testPhi, gamma=None, k=testK)
print("Solved for missing gamma_n:", missing_gamma_n)

# Case 2: Given values and solve for the missing r
missingK = functions.function21_59(phi=testPhi, gamma=testGamma, k=None)
print("Solved for missing K:", missingK)

# Case 3: Given values and solve for the missing phi_n
missingPhi = functions.function21_59(phi=None, gamma=testGamma, k=testK)
print("Solved for missing phi:", missingPhi)

