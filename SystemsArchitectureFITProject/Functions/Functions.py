from SystemsArchitectureFITProject.Functions.Functions_for_Functions import FunctionsFor21_85

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
            if n_1 is not None and P is not None:
                P = (n_1 * T * 1e6) / 77.6
                return P
            else:
                raise ValueError("You must provide n_1 and T to solve for P.")

        # If T is missing, calculate T (Temperature)
        elif T is None:
            if n_1 is not None and P is not None:
                T = (77.6 * P * 1e6)/ n_1
                return T
            else:
                raise ValueError("You must provide n_1 and P to solve for T")

        else:
            raise ValueError("One variable must be missing to solve for it.")

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
    