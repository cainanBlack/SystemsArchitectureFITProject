import numpy as np
from SystemsArchitectureFITProject.Functions.Functions import Functions

# Create an instance of the Functions class
functions = Functions()

# ============================================================================================== #
# Test Case 1: Solving for missing Hf
Hf = None
fx = 0.5
fy = 0.3
lambda_ = 0.6  # Wavelength in appropriate units
di = 2.0       # Propagation distance
Dx = 1.0       # Characteristic dimension in x
Dy = 1.5       # Characteristic dimension in y

# Solve for missing I_f
missing_Hf = functions.function21_32(Hf, fx, fy, lambda_, di, Dx, Dy)
print("Solved for missing H(f):", missing_Hf)

# ============================================================================================== #
# # # Test Case 2: Solving for missing lambda
fx = 0.5  # Spatial frequency in x-direction
fy = 0.3  # Spatial frequency in y-direction
di = 2.0  # Propagation distance
Dx = 1.0  # Characteristic dimension in x
Dy = 1.5  # Characteristic dimension in y
Hf = None
lambda_ = None


# # Solve for missing lambda
missing_lambda = functions.function21_32(fx=fx, fy=fy, di=di, Dx=Dx, Dy=Dy)
print("Solved for missing lambda:", missing_lambda)


# # # # ============================================================================================== #
# # # # # Test Case 3: Solving for missing di

Hf = 0.5  # Transfer function value
fx = 0.5       # Spatial frequency in x-direction
fy = 0.3       # Spatial frequency in y-direction
lambda_ = 0.6  # Wavelength in appropriate units
Dx = 1.0       # Characteristic dimension in x-direction
Dy = 1.5       # Characteristic dimension in y-direction

# # # Solve for missing Di
missing_di = functions.function21_32(Hf=Hf, fx=fx, fy=fy, lambda_=lambda_, Dx=Dx, Dy=Dy)
print("Solved for missing di:", missing_di)

# # # # ============================================================================================== #
# # # # # Test Case 4: Solving for missing Dx

Hf = 0.5       # Transfer function value
fx = 0.5       # Spatial frequency in x-direction
lambda_ = 0.6  # Wavelength in appropriate units
di = 2.0       # Propagation distance

# # # Solve for missing Dx
missing_Dx = functions.function21_32(Hf=Hf, fx=fx, lambda_=lambda_, di=di)
print("Solved for missing Dx:", missing_Dx)

# # # # ============================================================================================== #
# # # # # Test Case 5: Solving for missing Dy

Hf = 0.5       # Transfer function value
fy = 0.3       # Spatial frequency in y-direction
lambda_ = 0.6  # Wavelength in appropriate units
di = 2.0       # Propagation distance

# # Solve for missing lambda
missing_Dy = functions.function21_32(Hf=Hf, fy=fy, lambda_=lambda_, di=di)
print("Solved for missing Dy:", missing_Dy)

# # # ============================================================================================== #
# # # # Test Case 5: Solving for missing fx

# Hf = 0.5       # Transfer function value
# lambda_ = 0.6  # Wavelength in appropriate units
# di = 2.0       # Propagation distance
# Dx = 1.0       # Characteristic dimension in x-direction

# # # Solve for missing lambda
# missing_fx = functions.function21_32(Hf=Hf, lambda_=lambda_, di=di, Dx=Dx)
# print("Solved for missing fx:", missing_fx)

# # # # ============================================================================================== #
# # # # # Test Case 6: Solving for missing fy

# Hf = 0.5       # Transfer function value
# lambda_ = 0.6  # Wavelength in appropriate units
# di = 2.0       # Propagation distance
# Dy = 1.5       # Characteristic dimension in y-direction
# fx = None

# # # Solve for missing lambda
# missing_fy = functions.function21_32(Hf=Hf, lambda_=lambda_, di=di, Dx=Dx, Dy=Dy)
# print("Solved for missing fy:", missing_fy)

# # ============================================================================================== #
