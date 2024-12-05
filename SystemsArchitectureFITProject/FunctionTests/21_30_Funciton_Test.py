import numpy as np
from SystemsArchitectureFITProject.Functions.Functions import Functions

# Create an instance of the Functions class
functions = Functions()

# ============================================================================================== #
# Test Case 1: Solving for missing I_f
O_f = [complex(1, 1), complex(1.5, 2), complex(5, 5.5)]  # Example O(f)
H_f = [complex(1, 1), complex(2, 0), complex(1, -1)]     # Example H(f)

# Solve for missing I_f
missing_I_f = functions.function21_30(O_f=O_f, H_f=H_f)
print("Solved for missing I(f):", missing_I_f)

# ============================================================================================== #
# # Test Case 2: Solving for missing O_f
I_f = np.array([complex(1, 2), complex(3, 4), complex(5, 6)])   # Example values for I(f)
H_f = np.array([complex(1, 1), complex(2, 0), complex(1, -1)])  # Example values for H(f)

# Solve for missing O_f
missing_O_f = functions.function21_30(I_f=I_f, H_f=H_f)
print("Solved for missing O(f):", missing_O_f)

# ============================================================================================== #
# # Test Case 3: Solving for missing H_f
I_f = np.array([complex(1, 2), complex(3, 4), complex(5, 6)])       # Example I(f)
O_f = np.array([complex(1, 1), complex(1.5, 2), complex(5, 5.5)])   # Example O(f)

# Solve for missing H_f
missing_H_f = functions.function21_30(I_f=I_f, O_f=O_f)
print("Solved for missing H(f):", missing_H_f)

# ============================================================================================== #
