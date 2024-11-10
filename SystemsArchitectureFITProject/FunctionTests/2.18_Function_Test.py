from SystemsArchitectureFITProject.Functions.Functions import Functions

# Create an instance of the Functions class
functions = Functions()

# Case 1: Given values and solve for the missing d_i
f = 10  # Average frequency (f with bar)
x = 20  # Average value of x (x with bar)
lambda_ = 2  # Wavelength (known)
d_i = None  # d_i is missing

# Solve for the missing d_i using function2_18
missingDi = functions.function2_18(f=f, x=x, lambda_=lambda_)
print("Solved for missing d_i:", missingDi)

# Case 2: Given values and solve for the missing f
d_i = 5  # Distance (known)
f = None  # Frequency (f with bar) is missing

# Solve for the missing f using function2_18
missingF = functions.function2_18(x=x, lambda_=lambda_, d_i=d_i)
print("Solved for missing f:", missingF)

# Case 3: Given values and solve for the missing x
d_i = 5  # Distance (known)
lambda_ = 2  # Wavelength (known)
f = 10  # Average frequency (f with bar)

# Solve for the missing x using function2_18
missingX = functions.function2_18(f=f, lambda_=lambda_, d_i=d_i)
print("Solved for missing x:", missingX)

# Case 4: Given values and solve for the missing lambda
f = 10  # Frequency (f with bar)
x = 20  # Value of x (x with bar)
d_i = 5  # Distance (known)

# Solve for the missing lambda using function2_18
missingLambda = functions.function2_18(f=f, x=x, d_i=d_i)
print("Solved for missing lambda:", missingLambda)