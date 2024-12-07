import numpy as np
from SystemsArchitectureFITProject.Functions.Functions import Functions

# Create an instance of the Functions class
functions = Functions()

# ============================================================================================== #
# Test Case 1: Solving for missing Intensity

# Define parameters
Dx = 0.01  # Aperture width in x-direction (meters)
Dy = 0.01  # Aperture width in y-direction (meters)
lambda_ = 500e-9  # Wavelength (meters, e.g., 500 nm)
di = 1.0  # Distance (meters)
sxy = None

# Observation plane coordinates
x = np.linspace(-0.01, 0.01, 500)  # x-axis (meters)
y = np.linspace(-0.01, 0.01, 500)  # y-axis (meters)
X, Y = np.meshgrid(x, y)  # Create a 2D grid

# Calculate the intensity
intensity = functions.function21_33(sxy, X, Y, Dx, Dy, lambda_, di)
print("Solved for intensity:", intensity)

# ============================================================================================== #
# # # # Test Case 2: Solving for missing lambda
s = 1e4  # Example intensity value
x = 0.001  # x-coordinate (meters)
y = 0.001  # y-coordinate (meters)
Dx = 0.01  # Aperture width in x-direction (meters)
Dy = 0.01  # Aperture width in y-direction (meters)
di = 1.0  # Distance (meters)
lambda_ = None

# # # Solve for missing lambda
missing_lambda = functions.function21_33(sxy=s, x=x, y=y, Dx=Dx, Dy=Dy, lambda_=lambda_, di=di)
print("Solved for missing lambda:", missing_lambda)


# # # # ============================================================================================== #
# # # # Test Case 3: Solving for missing di
s = 1e4  # Example intensity value
x = 0.001  # x-coordinate (meters)
y = 0.001  # y-coordinate (meters)
Dx = 0.01  # Aperture width in x-direction (meters)
Dy = 0.01  # Aperture width in y-direction (meters)
lambda_ = 500e-9  # Wavelength (meters, e.g., 500 nm)
di = None  # Distance (meters)

# # # Solve for missing di
missing_di = functions.function21_33(sxy=s, x=x, y=y, Dx=Dx, Dy=Dy, lambda_=lambda_, di=di)
print("Solved for missing di:", missing_di)

# # # # ============================================================================================== #
# # # # Test Case 4: Solving for missing Dx
# Define known parameters
s = 1e4  # Example intensity value
x = 0.001  # x-coordinate (meters)
y = 0.001  # y-coordinate (meters)
Dy = 0.01  # Aperture height in y-direction (meters)
lambda_ = 500e-9  # Wavelength (meters, e.g., 500 nm)
di = 1.0  # Distance (meters)
Dx = None

missing_Dx = functions.function21_33(sxy=s, x=x, y=y, Dx=Dx, Dy=Dy, lambda_=lambda_, di=di)
print("Solved for missing Dx:", missing_Dx)

# # # # ============================================================================================== #
# # # # Test Case 5: Solving for missing Dy
# Define known parameters
s = 1e4  # Example intensity value
x = 0.001  # x-coordinate (meters)
y = 0.001  # y-coordinate (meters)
Dx = 0.01  # Aperture width in x-direction (meters)
lambda_ = 500e-9  # Wavelength (meters, e.g., 500 nm)
di = 1.0  # Distance (meters)
Dy = None

missing_Dy = functions.function21_33(sxy=s, x=x, y=y, Dx=Dx, Dy=Dy, lambda_=lambda_, di=di)
print("Solved for missing Dy:", missing_Dy)


# # # # ============================================================================================== #
# # # # Test Case 6: Solving for missing x
# Define known parameters
s = 1e4  # Example intensity value
y = 0.001  # y-coordinate (meters)
Dx = 0.01  # Aperture width in x-direction (meters)
Dy = 0.01  # Aperture height in y-direction (meters)
lambda_ = 500e-9  # Wavelength (meters, e.g., 500 nm)
di = 1.0  # Distance (meters)
x = None

missing_x = functions.function21_33(sxy=s, x=x, y=y, Dx=Dx, Dy=Dy, lambda_=lambda_, di=di)
print("Solved for missing x:", missing_x)


# # # # ============================================================================================== #
# # # # Test Case 6: Solving for missing y
# Define known parameters
s = 1e4  # Example intensity value
x = 0.001  # x-coordinate (meters)
Dx = 0.01  # Aperture width in x-direction (meters)
Dy = 0.01  # Aperture height in y-direction (meters)
lambda_ = 500e-9  # Wavelength (meters, e.g., 500 nm)
di = 1.0  # Distance (meters)
y = None

missing_y = functions.function21_33(sxy=s, x=x, y=y, Dx=Dx, Dy=Dy, lambda_=lambda_, di=di)
print("Solved for missing y:", missing_y)


