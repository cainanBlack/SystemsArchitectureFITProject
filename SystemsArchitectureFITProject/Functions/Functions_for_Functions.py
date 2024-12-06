from scipy.optimize import bisect
import math

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

class FunctionsFor21_19:
        
    @staticmethod
    def eulerFormula(theta_x):
        
        theta_radians = math.radians(theta_x)
        cos_theta = math.cos(theta_radians)
        sin_theta = math.sin(theta_radians)
        
        comp = complex(round(cos_theta, 4) , round(sin_theta, 4))
        return comp
        
        
