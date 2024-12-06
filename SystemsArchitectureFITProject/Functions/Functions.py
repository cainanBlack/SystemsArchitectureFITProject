from SystemsArchitectureFITProject.Functions.Functions_for_Functions import FunctionsFor21_30, FunctionsFor21_58, FunctionsFor21_59, FunctionsFor21_85
import math
import numpy as np

class Functions:
    
    # Solve for the missing variable in the equation f = x / (lambda * d_i), 
    # where f, x, lambda, and d_i may be known, and one is missing.
    # This function solves for the missing variable, either f, x, lambda, or d_i,
    # based on the provided known values. If any required variable is missing, 
    # the function will compute its value using the given equation.
    # @param f: Known average frequency (f with bar). If None, the function will solve for it.
    # @param x: Known average value of x (x with bar). If None, the function will solve for it.
    # @param lambda_: Known wavelength. If None, the function will solve for it.
    # @param d_i: Known distance or other parameter. If None, the function will solve for it.
    # @return: The computed value of the missing variable (f, x, lambda, or d_i).
    # @throws ValueError: If insufficient information is provided to solve for the missing variable.
    @staticmethod
    def function21_18(f=None, x=None, lambda_=None, d_i=None):
        
        # If f is missing, calculate f (average frequency)
        if f is None:
            if x is not None and lambda_ is not None and d_i is not None:
                f = x / (lambda_ * d_i)
                return f
            else:
                raise ValueError("You must provide x, lambda, and d_i to solve for f.")
        
        # If x is missing, calculate x (average value of x)
        elif x is None:
            if f is not None and lambda_ is not None and d_i is not None:
                x = f * lambda_ * d_i
                return x
            else:
                raise ValueError("You must provide f, lambda, and d_i to solve for x.")
        
        # If lambda is missing, calculate lambda
        elif lambda_ is None:
            if f is not None and x is not None and d_i is not None:
                lambda_ = x / (f * d_i)
                return lambda_
            else:
                raise ValueError("You must provide f, x, and d_i to solve for lambda.")
        
        # If d_i is missing, calculate d_i
        elif d_i is None:
            if f is not None and x is not None and lambda_ is not None:
                d_i = x / (f * lambda_)
                return d_i
            else:
                raise ValueError("You must provide f, x, and lambda to solve for d_i.")
        
        else:
            raise ValueError("One variable must be missing to solve for it.")

    # Solve for the missing variable in the equation I(f) = O(f) ∗ H(f).
    # where I(f), O(f), and H(f) may be known, and one is missing.
    # This function solves for the missing variable, either I(f) or  O(f) or H(f).
    # based on the provided known values. If any required variable is missing,
    # the function will compute its value using the given equation.
    # @param I_f: Known image spectrum. If None, the function will solve for it.
    # @param O_f: Known object radiant emittance. If None, the function will solve for it.
    # @param H_f: Known optical transfer function. If None, the function will solve for it.
    # @return: The computed value of the missing variable (I_f, O_f, H_f).
    # @throws ValueError: If insufficient information is provided to solve for the missing variable.
    @staticmethod
    def function21_30(I_f=None, O_f=None, H_f=None):

        # If  I_f is missing, calculate  I_f (image spectrum)
        if I_f is None:
            if O_f is not None and H_f is not None:
                return FunctionsFor21_30.compute_I_f(O_f, H_f)
            else:
                raise ValueError("You must provide O_f, and H_f to solve for I_f.")

        # If O_f is missing, calculate O_f (object radiant emittance)
        elif O_f is None:
            if I_f is not None and H_f is not None:
                return FunctionsFor21_30.compute_O_f(I_f, H_f)
            else:
                raise ValueError("You must provide I_f and H_f to solve for O_f.")

        # If H_f is missing, calculate H_f (optical transfer function)
        elif H_f is None:
            if I_f is not None and O_f is not None:
                return FunctionsFor21_30.compute_H_f(I_f, O_f)
            else:
                raise ValueError("You must provide I_f and O_f to solve for H_f")

        else:
            raise ValueError("One variable must be missing to solve for it.")

    # Solve for the missing variable in the equation n_1 =  (77.6P / T) * 10^(-6)
    # where n_1, P, and T may be known, and one is missing.
    # This function solves for the missing variable, either n_1, P, or T,
    # based on the provided known values. If any required variable is missing,
    # the function will compute its value using the given equation.
    # @param n_1: Known index of refraction. If None, the function will solve for it.
    # @param P: Known function of pressure (P, in millibars). If None, the function will solve for it.
    # @param T: Known funtion of temperature (T, in Kelvin). If None, the function will solve for it.
    # @return: The computed value of the missing variable (n_1, P, T).
    # @throws ValueError: If insufficient information is provided to solve for the missing variable.
    @staticmethod
    def function21_44(n_1=None, P=None, T=None):

        # If n_1 is missing, calculate n_1 (index of refraction)
        if n_1 is None:
            if P is not None and T is not None:
                n_1 = ((77.6 * P) / T) * 1e-6
                return n_1
            else:
                raise ValueError("You must provide P, and T to solve for n_1.")

        # If P is missing, calculate P (Pressure)
        elif P is None:
            if n_1 is not None and T is not None:
                P = (n_1 * T * 1e6) / 77.6
                return P
            else:
                raise ValueError("You must provide n_1 and T to solve for P.")

        # If T is missing, calculate T (Temperature)
        elif T is None:
            if n_1 is not None and P is not None:
                T = (77.6 * P) / (n_1 * 1e6)
                return T
            else:
                raise ValueError("You must provide n_1 and P to solve for T")

        else:
            raise ValueError("One variable must be missing to solve for it.")
        
    # Solves for the missing variable: phi_n, r, or gamma_n.
    # Equation Components:
    #     𝛤𝑛(𝑟) = ∫ 𝑑𝑘⃗ 𝛷𝑛(𝑘⃗)𝑒(−𝑗𝑘⃗ ∗𝑟)
    # @param phi_n 𝛷𝑛(𝑘⃗)
    # @param r The position vector in real space
    # @param gamma_n The value of Γₙ(r)
    # @return: The computed value of the equation
    @staticmethod
    def function21_58(phi_n=None, r=None, gamma_n=None):

        if gamma_n is None:
            if phi_n is not None and r is not None:
                return FunctionsFor21_58.solveForGamma_N(phi_n, r)
            else:
                raise ValueError("phi_n and r are required to solve for gamma_n.")

        elif phi_n is None:
            if r is not None and gamma_n is not None:
                return FunctionsFor21_58.solveForPhi_n(gamma_n, r)
            else:
                raise ValueError("r and gamma_n are required to solve for phi_n.")

        elif r is None:
            if phi_n is not None and gamma_n is not None:
                return FunctionsFor21_58.solveForR(phi_n, gamma_n)
            else:
                raise ValueError("phi_n and gamma_n are required to solve for r.")
            
        else:
            raise ValueError("Invalid missing variable. Choose 'gamma_n', 'phi_n', or 'r'.")

    # Solves for and missing value in 𝛷𝑛(𝑘⃗ ) = 1\(2𝜋)3 ∫ 𝑑𝑟𝛤𝑛(𝑟)𝑒(𝑗𝑘⃗ ∗𝑟)n
    # @param phi Target value for Φ_n(k) (set to None if solving for Φ_n(k)).
    # @param gamma Grid function for Γ_n(r) (set to None if solving for Γ_n(r)).
    # @param k Array [kx, ky, kz] (set to None if solving for k).
    # @param grid_points Number of points along each dimension of the grid.
    # @param bounds Tuple (min, max) defining the spatial and Fourier bounds.
    # @return The solved variable (Φ_n(k), Γ_n(r), or k).
    @staticmethod
    def function21_59(phi=None, gamma=None, k=None, grid_points=32, bounds=(-10, 10)):
        if phi is None:
            # Solve for Φ_n(k)
            if k is None:
                raise ValueError("Missing k for solving Φ_n(k).")
            gamma_grid = FunctionsFor21_59.gamma_n_grid(grid_points, bounds)
            phi_k = FunctionsFor21_59.compute_phi_fft(gamma_grid, bounds)
            kx, ky, kz = k
            k_values = np.linspace(bounds[0], bounds[1], grid_points)
            i = (np.abs(k_values - kx)).argmin()
            j = (np.abs(k_values - ky)).argmin()
            k = (np.abs(k_values - kz)).argmin()
            return phi_k[i, j, k]

        elif gamma is None:
            # Solve for Γ_n(r)
            if k is None:
                raise ValueError("Missing k for solving Γ_n(r).")
            phi_grid = FunctionsFor21_59.compute_phi_fft(FunctionsFor21_59.gamma_n_grid(grid_points, bounds), bounds)
            gamma_grid = FunctionsFor21_59.compute_gamma_fft(phi_grid, bounds)
            r_values = np.linspace(bounds[0], bounds[1], grid_points)
            x, y, z = k
            i = (np.abs(r_values - x)).argmin()
            j = (np.abs(r_values - y)).argmin()
            k = (np.abs(r_values - z)).argmin()
            return gamma_grid[i, j, k]

        elif k is None:
            # Solve for k
            gamma_grid = FunctionsFor21_59.gamma_n_grid(grid_points, bounds)
            return FunctionsFor21_59.solve_for_k(phi, gamma_grid, bounds)

        else:
            raise ValueError("One variable must be None to solve for it.")

    # Solves for any missing variable: phi, a_i, rho, or theta.
    # Depending on which variable is missing, it either computes the missing value 
    # or raises a ValueError if the required variables are not provided.
    # @param phi: The value of phi. If None, the function will compute it.
    # @param a: A list of coefficients a_i. If any coefficient is None, the function will solve for the missing coefficient.
    # @param rho: The value of rho. If None, the function will solve for the missing rho.
    # @param theta: The value of theta. If None, the function will solve for the missing theta.
    # @param Z: A list of lambda functions used to compute the Z_i values for the variables. Each Z_i is a function of rho and theta.
    # @return: The computed value of the missing variable (phi, a_i, rho, or theta).
    #          If all variables are provided, the value of phi is returned.
    # @throws ValueError: If insufficient information is provided to solve for a missing variable.
    @staticmethod
    def function21_85(phi=None, a=None, rho=None, theta=None, Z=None):
        """
        Solves for any missing variable: phi, a_i, rho, or theta.
        """

        # Handle missing `phi`
        if phi is None:
            if a is None or rho is None or theta is None or Z is None:
                raise ValueError("Need a, rho, theta, and Z to compute phi")
            return FunctionsFor21_85.computePhi(rho, theta, a, Z)

        # Handle missing `a_i`
        if a is not None and any(a_i is None for a_i in a):
            if phi is None or rho is None or theta is None or Z is None:
                raise ValueError("Need phi, rho, theta, and Z to solve for a_i")
            return FunctionsFor21_85.solveForA(phi, rho, theta, Z, a)

        # Handle missing `rho`
        if rho is None:
            if phi is None or a is None or theta is None or Z is None:
                raise ValueError("Need phi, a, theta, and Z to solve for rho")
            return FunctionsFor21_85.solveForRho(phi, a, theta, Z)
        
        # Handle missing `theta`
        if theta is None:
            if phi is None or a is None or rho is None or Z is None:
                raise ValueError("Need phi, a, rho, and Z to solve for theta")
            return FunctionsFor21_85.solveForTheta(phi, a, rho, Z)

        # If everything is provided, return phi
        return phi
    
    # Solves for the missing variable in the SNRD equation.
    # Equation components:
    #     SNRD(f) = K_bar * |H(f)| * |O_n(f)| / sqrt(K_bar + (K_bar)^2 * |O_n(f)|^2) * sqrt(var{H(f)} + P * sigma_n^2)    
    # @param f Frequency value (not used directly in calculation)
    # @param h_f Signal response at frequency f (optional, None if unknown)
    # @param Noise signal at frequency f (optional, None if unknown)
    # @param Scaling factor (optional, None if unknown)
    # @param Variance of h_f (optional, None if unknown)
    # @param p_sigma_n_squared Noise power (optional, None if unknown)
    # @return The calculated value of the missing variable.
    @staticmethod
    def function21_104(h_f=None, o_n_f=None, k_bar=None, var_h_f=None, p_sigma_n_squared=None):
        
        # If h_f is missing, solve for it
        if h_f is None:
            if all(v is not None for v in [o_n_f, k_bar, var_h_f, p_sigma_n_squared]):
                numerator = k_bar * abs(o_n_f)
                denominator = math.sqrt(k_bar + k_bar**2 * abs(o_n_f)**2) * math.sqrt(var_h_f + p_sigma_n_squared)
                return numerator / denominator
            else:
                raise ValueError("Insufficient data to solve for h_f.")
        
        # If o_n_f is missing, solve for it
        elif o_n_f is None:
            if all(v is not None for v in [h_f, k_bar, var_h_f, p_sigma_n_squared]):
                numerator = k_bar * abs(h_f)
                denominator = math.sqrt(k_bar + k_bar**2 * (10**-10)) * math.sqrt(var_h_f + p_sigma_n_squared)
                result = numerator / denominator
                
                # Now solve for o_n(f)
                o_n_f_result = abs(result)
                return o_n_f_result
            else:
                raise ValueError("Insufficient data to solve for o_n_f.")
        
        # If k_bar is missing, solve for it
        elif k_bar is None:
            if all(v is not None for v in [h_f, o_n_f, var_h_f, p_sigma_n_squared]):
                numerator = abs(h_f) * abs(o_n_f)
                denominator = math.sqrt(var_h_f + p_sigma_n_squared)
                return numerator / denominator
            else:
                raise ValueError("Insufficient data to solve for k_bar.")
        
        # If var_h_f is missing, solve for it
        elif var_h_f is None:
            if all(v is not None for v in [h_f, o_n_f, k_bar, p_sigma_n_squared]):
                numerator = abs(h_f) * abs(o_n_f)
                denominator = math.sqrt(k_bar + k_bar**2 * abs(o_n_f)**2) * math.sqrt(p_sigma_n_squared)
                return numerator / denominator
            else:
                raise ValueError("Insufficient data to solve for var_h_f.")
        
        # If p_sigma_n_squared is missing, solve for it
        elif p_sigma_n_squared is None:
            if all(v is not None for v in [h_f, o_n_f, k_bar, var_h_f]):
                numerator = abs(h_f) * abs(o_n_f)
                denominator = math.sqrt(k_bar + k_bar**2 * abs(o_n_f)**2) * math.sqrt(var_h_f)
                return numerator / denominator
            else:
                raise ValueError("Insufficient data to solve for p_sigma_n_squared.")
        
        else:
            raise ValueError("No missing variable to solve for.")
        