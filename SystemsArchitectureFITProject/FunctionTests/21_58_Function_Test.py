from SystemsArchitectureFITProject.Functions.Functions import Functions
import numpy as np

# Create an instance of the Functions class
functions = Functions()

testGamma_n = 5 + 2j  # Γₙ(r)
testR = [1.0, 0.5]    # Position vector r
testK = (1,1,1)

# Case 1: Given values and solve for the missing gamma_n
missing_gamma_n = functions.function21_58(r=testR, k=testK)
print("Solved for missing n_1:", missing_gamma_n)

# Case 2: Given values and solve for the missing r
missingR = functions.function21_58(gamma_n=testGamma_n, k=testK)
print("Solved for missing R:", missingR)

# Case 3: Given values and solve for the missing k
missingPhi_n = functions.function21_58(r=testR, gamma_n=testGamma_n)
print("Solved for missing phi_n:", missingPhi_n)

