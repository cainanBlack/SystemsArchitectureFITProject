from SystemsArchitectureFITProject.Functions.Functions import Functions

# Create an instance of the Functions class
functions = Functions()

# Solve for h
test_h = 500
test_v = 30
test_A = -53
test_Cn2_value = 1e-14
h_missing = Functions.function21_52(h=None, v=test_v, A=test_A, Cn2_value=test_Cn2_value)
print("Solved h = ", h_missing)

# Solve for v
v_missing = Functions.function21_52(h=test_h, v=None, A=test_A, Cn2_value=test_Cn2_value)
print("Solved v = ", v_missing)

# Solve for A
A_missing = Functions.function21_52(h=test_h, v=test_v, A=None, Cn2_value=test_Cn2_value)
print("Solved A = ", A_missing)

# Solve for A
A_missing = Functions.function21_52(h=test_h, v=test_v, A=test_A, Cn2_value=None)
print("Solved Cn2 = ", A_missing)