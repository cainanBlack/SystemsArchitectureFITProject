from scipy.optimize import bisect
import math
import numpy as np
from scipy.integrate import nquad
from scipy.fft import fftn, ifftn, fftshift
from scipy.optimize import minimize

class FunctionsFor21_19:
        
    @staticmethod
    def eulerFormula(theta_x):
        
        theta_radians = math.radians(theta_x)
        cos_theta = math.cos(theta_radians)
        sin_theta = math.sin(theta_radians)
        
        comp = complex(round(cos_theta, 4) , round(sin_theta, 4))
        return comp
    
class FunctionsFor21_30:
    @staticmethod
    def compute_I_f(O_f, H_f):
        """
        Solves for I(f) in the equation I(f) = O(f) * H(f).

        Parameters:
            O_f (numpy.ndarray): The original signal in the frequency domain.
            H_f (numpy.ndarray): The transfer function in the frequency domain.

        Returns:
            numpy.ndarray: The output signal I(f) in the frequency domain.
        """
        O_f = np.array(O_f)
        H_f = np.array(H_f)

        # Solve for I(f) using element-wise multiplication
        I_f = O_f * H_f

        return I_f

    @staticmethod
    def compute_O_f(I_f, H_f):
        """
        Solves for O(f) in the equation I(f) = O(f) * H(f).

        Parameters:
            I_f (numpy.ndarray): The output signal in the frequency domain.
            H_f (numpy.ndarray): The transfer function in the frequency domain.

        Returns:
            numpy.ndarray: The original signal O(f) in the frequency domain.
        """

        I_f = np.array(I_f)
        H_f = np.array(H_f)

        # Ensure H_f does not have zero values to avoid division by zero
        H_f_safe = np.where(H_f == 0, np.nan, H_f) # Replace zeros with NaN

        # Solve for O(f)
        O_f = I_f / H_f_safe
        return O_f

    @staticmethod
    def compute_H_f(I_f, O_f):
        """
        Solves for H(f) in the equation I(f) = O(f) * H(f).

        Parameters:
            I_f (numpy.ndarray): The output signal in the frequency domain.
            O_f (numpy.ndarray): The original signal in the frequency domain.

        Returns:
            numpy.ndarray: The transfer function H(f) in the frequency domain.
        """

        I_f = np.array(I_f)
        O_f = np.array(O_f)

        # Ensure O_f does not have zero values to avoid division by zero
        O_f_safe = np.where(O_f == 0, np.nan, O_f)  # Replace zeros with NaN

        # Solve for H(f)
        H_f = I_f / O_f_safe

        return H_f

class FunctionsFor21_85:
    
    # Compute the value of phi based on the formula:
    # phi(rho, theta) = sum(a_i * Z_i(rho, theta)) for each valid a_i.
    # This function computes the value of phi using the given `a` coefficients and `Z_i` functions. It skips any `None` values in `a`.
    # @param rho: The value of rho to be used in the computation.
    # @param theta: The value of theta to be used in the computation.
    # @param a: A list of coefficients a_i used in the computation. The function skips any None values.
    # @param Z: A list of lambda functions representing Z_i(rho, theta) to be evaluated for each coefficient.
    # @return: The computed value of phi.
    @staticmethod
    def computePhi(rho, theta, a, Z):
        """
        Compute the value of phi based on the formula:
        phi(Rho, Theta) = sum(a_i * Z_i(rho, theta))
        Skip any `None` values in `a`.
        """
        return sum(a_i * Z_i(rho, theta) for a_i, Z_i in zip(a, Z) if a_i is not None)

    # Solve for the missing rho using the bisection method.
    # This method finds the value of rho that makes the computed value of phi equal to the given phi value.
    # @param phi: The value of phi that we want to match by solving for rho.
    # @param a: A list of coefficients a_i.
    # @param theta: The value of theta to be used in the computation.
    # @param Z: A list of lambda functions representing Z_i(rho, theta).
    # @return: The computed value of rho that satisfies the equation phi = sum(a_i * Z_i(rho, theta)).
    # @throws ValueError: If phi, a, theta, or Z is None.
    @staticmethod
    def solveForRho(phi, a, theta, Z):
        """
        Solve for the missing rho using the bisection method.
        """
        def func(rho_guess):
            return FunctionsFor21_85.computePhi(rho_guess, theta, a, Z) - phi
        
        # Bisection method to find root of func(rho) = 0
        return bisect(func, 0.001, 100.0, xtol=1e-6)
    
    # Solve for the missing theta using the bisection method.
    # This method finds the value of theta that makes the computed value of phi equal to the given phi value.
    # @param phi: The value of phi that we want to match by solving for theta.
    # @param a: A list of coefficients a_i.
    # @param rho: The value of rho to be used in the computation.
    # @param Z: A list of lambda functions representing Z_i(rho, theta).
    # @return: The computed value of theta that satisfies the equation phi = sum(a_i * Z_i(rho, theta)).
    # @throws ValueError: If phi, a, rho, or Z is None.
    @staticmethod
    def solveForTheta(phi, a, rho, Z):
        """
        Solve for the missing theta using the bisection method.
        """
        def func(thetaGuess):
            return FunctionsFor21_85.computePhi(rho, thetaGuess, a, Z) - phi
        
        # Bisection method to find root of func(theta) = 0
        return bisect(func, 0.001, 100.0, xtol=1e-6)
    
    # Solve for missing a_i coefficients in the equation:
    # phi = sum(a_i * Z_i(rho, theta)) for all i.
    # This function solves for each missing a_i by isolating each one in the equation and solving for it.
    # @param phi: The value of phi that we want to match by solving for the coefficients a_i.
    # @param rho: The value of rho to be used in the computation.
    # @param theta: The value of theta to be used in the computation.
    # @param Z: A list of lambda functions representing Z_i(rho, theta).
    # @param a: A list of coefficients a_i, where missing values (None) will be solved for.
    # @return: The list of coefficients a_i with the missing values filled in.
    # @throws ValueError: If the number of elements in a and Z do not match.
    @staticmethod
    def solveForA(phi, rho, theta, Z, a):
        """
        Solve for missing a_i coefficients in the equation:
        phi = sum(a_i * Z_i(rho, theta)) for all i
        """
        # Check if the lengths of a and Z match
        if len(a) != len(Z):
            raise ValueError("The number of elements in 'a' and 'Z' must be the same.")
        
        # Solve for each missing a_i
        for i, a_i in enumerate(a):
            if a_i is None:
                # Solve for this a_i using the equation:
                # a_i = (phi - sum(a_j * Z_j(rho, theta) for all j != i)) / Z_i(rho, theta)
                a[i] = (phi - sum(a_j * Z[j](rho, theta) for j, a_j in enumerate(a) if j != i)) / Z[i](rho, theta)
        
        return a        
        
class FunctionsFor21_58:

    # Solves for the missing value of gamma_n.
    # @param phi_n 𝛷𝑛(𝑘⃗)
    # @param r The position vector in real space
    # @return: The computed value of gamma_n
    @staticmethod
    def solveForGamma_N(phi_n, r):
            def integrand(*k):
                k = np.array(k)
                return phi_n(k) * np.exp(-1j * np.dot(k, r))

            # Define integration limits (assuming infinite)
            ndim = len(r)
            limits = [(-np.inf, np.inf)] * ndim

            # Perform the integration
            enter = integrand
            result = nquad(enter, limits)
            return result
    
    # Solves for the missing value of phi_n.
    # @param r The position vector in real space
    # @param gamma_n The value of Γₙ(r)
    # @return: The computed value of phi_n
    @staticmethod
    def solveForPhi_n(gamma_n, r):
            def phi_n(k):
                k = np.array(k)
                return gamma_n * np.exp(1j * np.dot(k, r))
            result = phi_n([0.2, 0.1])
            return result
     
    # Solves for the missing value of r.
    # @param phi_n 𝛷𝑛(𝑘⃗)
    # @param r The position vector in real space
    # @return: The computed value of r
    @staticmethod
    def solveForR(phi_n, gamma_n):
        def error_function(r_guess):
            r_guess = np.array(r_guess)

            def integrand(*k):
                k = np.array(k)
                return phi_n(k) * np.exp(-1j * np.dot(k, r_guess))

            # Integration limits (assuming infinite)
            ndim = len(r_guess)
            limits = [(-np.inf, np.inf)] * ndim
            result, _ = nquad(integrand, limits)
            return np.abs(result - gamma_n)

        # Initial guess for r
        ndim = len(gamma_n) if isinstance(gamma_n, (list, np.ndarray)) else 2
        initial_guess = np.zeros(ndim)

        result = minimize(error_function, initial_guess)
        return result.x

# Solves for the missing value of r.
# @param phi_n 𝛷𝑛(𝑘⃗)
# @param r The position vector in real space
# @return: The computed value of r    
class FunctionsFor21_59:
    
    @staticmethod
    def gamma_n_grid(grid_points, bounds):
        x = np.linspace(bounds[0], bounds[1], grid_points)
        y = np.linspace(bounds[0], bounds[1], grid_points)
        z = np.linspace(bounds[0], bounds[1], grid_points)
        xx, yy, zz = np.meshgrid(x, y, z, indexing="ij")
        r_squared = xx**2 + yy**2 + zz**2
        return np.exp(-r_squared)  # Example: Gaussian
    
    @staticmethod
    def compute_phi_fft(gamma_grid, bounds):
        """
        Compute Φ_n(k) using FFT.
        """
        dV = (bounds[1] - bounds[0]) / gamma_grid.shape[0]  # Volume element
        phi_k = fftshift(fftn(gamma_grid)) * dV
        return phi_k

    @staticmethod
    def compute_gamma_fft(phi_grid, bounds):
        """
        Compute Γ_n(r) using the inverse FFT.
        """
        dV = (bounds[1] - bounds[0]) / phi_grid.shape[0]  # Volume element
        gamma_r = ifftn(fftshift(phi_grid)) * phi_grid.size * dV
        return np.real(gamma_r)

    @staticmethod
    def solve_for_k(phi_target, gamma_grid, bounds):
        """
        Solve for k numerically using optimization.
        """
        phi_k = FunctionsFor21_59.compute_phi_fft(gamma_grid, bounds)
        k_values = np.linspace(bounds[0], bounds[1], gamma_grid.shape[0])

        def error_function(k):
            kx, ky, kz = k
            i = (np.abs(k_values - kx)).argmin()
            j = (np.abs(k_values - ky)).argmin()
            k = (np.abs(k_values - kz)).argmin()
            return np.abs(phi_k[i, j, k] - phi_target)

        initial_guess = [0, 0, 0]
        result = minimize(
            error_function,
            initial_guess,
            bounds=[(bounds[0], bounds[1])] * 3,
        )
        return result.x
