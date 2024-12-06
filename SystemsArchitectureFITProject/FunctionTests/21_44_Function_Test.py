from SystemsArchitectureFITProject.Functions.Functions import Functions

# Create an instance of the Functions class
functions = Functions()

# Case 1: Given values and solve for the missing n_1
P = 1007.33 # pressure (P, in millibars)
T = 300     # temperature (T, in Kelvin).
n_1 = None  # n_1 is missing

# Solve for the missing n_1 using function21_44
missing_n1 = functions.function21_44(T=T, P=P)
print("Solved for missing n_1:", missing_n1)

# Case 2: Given values and solve for the missing P
n_1 = 0.0003 # index of Refractivity (known)
T = 300      # temperature (T, in Kelvin).
P = None     # Pressure (P) is missing

# Solve for the missing P using function21_44
missingP = functions.function21_44(n_1=n_1, T=T)
print("Solved for missing P:", missingP)

# Case 3: Given values and solve for the missing T
n_1 = 0.0003  # index of Refractivity (known)
P = 1007.33   # pressure (known)
T =  None     # Temperature (T) is missing

# Solve for the missing T using function21_44
missingT = functions.function21_44(n_1=n_1, P=P)
print("Solved for missing T:", missingT)