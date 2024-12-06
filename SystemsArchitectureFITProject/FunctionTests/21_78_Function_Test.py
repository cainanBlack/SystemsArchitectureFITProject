from SystemsArchitectureFITProject.Functions.Functions import Functions

# Create an instance of the Functions class
functions = Functions()

# Solve for I_SE (scalar inputs)
I_SE = functions.function21_78(O=2, H_SE=1, H_0=2, I_SE=None)
print("Result I_SE:", I_SE)  # Output: 4

# Solve for O
O = functions.function21_78(I_SE=12, H_SE=1, H_0=3, O=None)
print("Result O:", O)  # Output: 4

# Solve for H_SE
H_SE = functions.function21_78(I_SE=12, O=4, H_0=3, H_SE=None)
print("Result H_SE:", H_SE)  # Output: 1

# Solve for H_0
H_0 = functions.function21_78(I_SE=12, O=4, H_SE=1, H_0=None)
print("Result H_0:", H_0)  # Output: 3
