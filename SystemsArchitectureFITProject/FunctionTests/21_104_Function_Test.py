from SystemsArchitectureFITProject.Functions.Functions import Functions

functions = Functions()

# Test for missing h_f
print("Test: Solve for h_f")
o_n_f = 2  # Noise signal at frequency f
k_bar = 5  # Scaling factor
var_h_f = 4  # Variance of H(f)
p_sigma_n_squared = 5  # Noise power
h_f_result = functions.function21_104(o_n_f=o_n_f, k_bar=k_bar, var_h_f=var_h_f, p_sigma_n_squared=p_sigma_n_squared)
print(f"Calculated h_f: {h_f_result}\n")

# Test for missing o_n_f
print("Test: Solve for o_n_f")
h_f = 10  # Signal response at frequency f
k_bar = 5  # Scaling factor
var_h_f = 4  # Variance of H(f)
p_sigma_n_squared = 5  # Noise power
o_n_f_result = functions.function21_104(h_f=h_f, k_bar=k_bar, var_h_f=var_h_f, p_sigma_n_squared=p_sigma_n_squared)
print(f"Calculated o_n_f: {o_n_f_result}\n")

# Test for missing k_bar
print("Test: Solve for k_bar")
h_f = 10  # Signal response at frequency f
o_n_f = 2  # Noise signal at frequency f
var_h_f = 4  # Variance of H(f)
p_sigma_n_squared = 5  # Noise power
k_bar_result = functions.function21_104(h_f=h_f, o_n_f=o_n_f, var_h_f=var_h_f, p_sigma_n_squared=p_sigma_n_squared)
print(f"Calculated k_bar: {k_bar_result}\n")

# Test for missing var_h_f
print("Test: Solve for var_h_f")
h_f = 10  # Signal response at frequency f
o_n_f = 2  # Noise signal at frequency f
k_bar = 5  # Scaling factor
p_sigma_n_squared = 5  # Noise power
var_h_f_result = functions.function21_104(h_f=h_f, o_n_f=o_n_f, k_bar=k_bar, p_sigma_n_squared=p_sigma_n_squared)
print(f"Calculated var_h_f: {var_h_f_result}\n")

# Test for missing p_sigma_n_squared
print("Test: Solve for p_sigma_n_squared")
h_f = 10  # Signal response at frequency f
o_n_f = 2  # Noise signal at frequency f
k_bar = 5  # Scaling factor
var_h_f = 4  # Variance of H(f)
p_sigma_n_squared_result = functions.function21_104(h_f=h_f, o_n_f=o_n_f, k_bar=k_bar, var_h_f=var_h_f)
print(f"Calculated p_sigma_n_squared: {p_sigma_n_squared_result}\n")