from SystemsArchitectureFITProject.Functions.Functions import Functions
import numpy as np

# Create an instance of the Functions class
functions = Functions()

# Values for Test
def testPhi_n(k):
    # Example Φₙ(k⃗): Gaussian form
    return np.exp(-np.linalg.norm(k)**2)

testGamma_n = 5 + 2j  # Γₙ(r)
testR = [1.0, 0.5]    # Position vector r

# Case 1: Given values and solve for the missing gamma_n
missing_gamma_n = functions.function21_58(phi_n=testPhi_n, r=testR)
print("Solved for missing n_1:", missing_gamma_n)

# Case 2: Given values and solve for the missing r
missingR = functions.function21_58(phi_n=testPhi_n, gamma_n=testGamma_n)
print("Solved for missing R:", missingR)

# Case 3: Given values and solve for the missing phi_n
missingPhi_n = functions.function21_58(r=testR, gamma_n=testGamma_n)
print("Solved for missing phi_n:", missingPhi_n)

