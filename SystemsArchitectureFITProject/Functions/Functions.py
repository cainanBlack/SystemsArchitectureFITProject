from SystemsArchitectureFITProject.Functions.Functions_for_Functions import FunctionsFor21_19, FunctionsFor21_85
import math

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

    # Solve for the missing variable in the equation gpf = wp * e^(j(theta_x)) 
    # @param gpf: Known generalized pupil function. If None, the function will solve for it.
    # @param wp: Known amplitude function. If None, the function will solve for it.
    # @param theta_x: Known phase term. If None, the function will solve for it.
    # @return: The computed value of the missing variable
    # @throws ValueError: If insufficient information is provided to solve for the missing variable.    
    @staticmethod
    def function21_19(gpf=None, wp=None, theta_x=None):
        
        #if gpf is missing, solve for gpf
        if gpf is None:
            if wp is not None and theta_x is not None:
                gpf = wp * FunctionsFor21_19.eulerFormula(theta_x)
                return gpf
            else:
                raise ValueError("You must provide wp and theta to solve for gpf")
        
        #if wp is missing, solve for wp
        if wp is None:
            if gpf is not None and theta_x is not None:
                wp = gpf / FunctionsFor21_19.eulerFormula(theta_x)
                return wp
            else:
                raise ValueError("You must provide gpf and theta to solve for wp")
            
        #if theta_x is missing, solve for theta_x
        if theta_x is None:
            if gpf is not None and wp is not None:
                theta_x = (math.log(gpf/wp) / complex(0, 1))
                return theta_x
            else:
                raise ValueError("You must provide gpf and wp to solve for theta_x")   
            
    # Solve for the missing variable or coordinates in the equations delta_x = (lambda_ * z) / d_x and delta_y = (lambda_ * z) / d_y
    # Where delta_x, delta_y, lambda, z, d_y, and d_x may be known and one is missing. delta_x and delta_y may both be missing if the rest
    # of the variables are known. This function solves for each variable based on the provided known values. If any required variable(s) is missing,
    # the funciton will compute its value using the given equation.
    # @param delta_x: the x coordinate of the given rectangular aperture. If None, the function will solve for it.
    # @param delta_y: the y coordinate of the given rectangular aperture. If None, the function will solve for it.
    # @param lambda_: Known wavelength. If None, the function will solve for it.
    # @param z: Known separation distance between object plane and entrance pupil plane of the imaging system. If None, the function will solve for it.
    # @param d_x: Known length of the x direction aperture width. If None, the function will solve for it.
    # @param d_y: Known length of the y direction aperture width. If None, the function will solve for it.
    # @return: The computed value of the missing variable(s)
    # @throws ValueError: If insufficient information is provided to solve for the missing variable.
    @staticmethod 
    def function21_38(delta_x=None, delta_y=None, lambda_=None, z=None, d_x=None, d_y=None):
        
        # If both coordinates are missing, calculate them 
        if delta_x is None and delta_y is None:
            if lambda_ is not None and z is not None and d_x is not None and d_y is not None:
                delta_x = (lambda_ * z) / d_x
                delta_y = (lambda_ * z) / d_y
                return delta_x, delta_y
            else:
                raise ValueError("You must provide lambda, z, and d_x to solve for delta_x and lambda, z, and d_y to solve for delta_y")
        
        # If delta_x is missing and delta_y is known, calculate delta_x
        elif delta_x is None and delta_y is not None and d_x is not None and z is not None:
            if lambda_ is not None:
                delta_x = (lambda_ * z) / d_x
                return delta_x, delta_y
            else:
                raise ValueError("You must provide lambda, z, and d_x to solve for delta_x")
        
        # If delta_y is missing and delta_x is known, calculate delta_y
        elif delta_y is None and delta_x is not None and d_y is not None and z is not None:
            if lambda_ is not None:
                delta_y = (lambda_ * z) / d_y
                return delta_x, delta_y
            else:
                raise ValueError("You must provide lambda, z, and d_y to solve for delta_y")
        
        # If lambda_ is missing and delta_x is known, calculate lambda_
        elif  lambda_ is None and delta_x is not None:
            if z is not None and d_x is not None:
                lambda_ = (delta_x * d_x) / z
                return lambda_
            else:
                raise ValueError("You must provide delta_x, d_x, and z to solve for lambda")
        
        # If lambda_ is missing and delta_y is known, calculate lambda_    
        elif lambda_ is None and delta_y is not None:
            if z is not None and d_y is not None:
                lambda_ = (delta_y * d_y) / z
                return lambda_
            else:
                raise ValueError("You must provide delta_y, d_y, and z to solve for lambda")
        
        # If z is missing and delta_x is known, calculate z    
        elif z is None and delta_x is not None:
            if lambda_ is not None and d_x is not None:
                z = (delta_x * d_x) / lambda_
                return z
            else:
                raise ValueError("You must provide delta_x, d_x, and lambda to solve for z")
        
        # If z is missing and delta_y is known, calculate z      
        elif z is None and delta_y is not None:
            if lambda_ is not None and d_y is not None:
                z = (delta_y * d_y) / lambda_
                return z
            else:
                raise ValueError("You must provide delta_y, d_y, and lambda to solve for z")
        
        # If d_x is missing, calculate d_x
        elif d_x is None and delta_x is not None:
            if z is not None and lambda_ is not None:
                d_x = (lambda_ * z) / delta_x
                return d_x
            else:
                raise ValueError("You must provide lambda, z, and delta_x to solve for d_x")
        
        # If d_y is missing, calculate d_y
        elif d_y is None and delta_y is not None:
            if z is not None and lambda_ is not None:
                d_y = (lambda_ * z) / delta_y
                return d_y
            else:
                raise ValueError("You must provide lambda, z, and delta_y to solve for d_y")
        
        else:
            raise ValueError("Invalid combination of missing variables")
    
    # Solve for the missing variable or coordinates in the equations delta_theta_d_x = lambda_ / d_x and delta_theta_d_y = lambda_ / d_y
    # Where delta_theta_d_x, delta_theta_d_y, lambda, d_y, and d_x may be known and one is missing. delta_theta_d_x and delta_theta_d_y may both be missing if the rest
    # of the variables are known. This function solves for each variable based on the provided known values. If any required variable(s) is missing,
    # the funciton will compute its value using the given equation.
    # @param delta_theta_d_x: the x coordinate of the given rectangular aperture based on angular properties. If None, the function will solve for it.
    # @param delta_theta_d_y: the y coordinate of the given rectangular aperture based on angular properties. If None, the function will solve for it.
    # @param lambda_: Known wavelength. If None, the function will solve for it.
    # @param d_x: Known length of the x direction aperture width. If None, the function will solve for it.
    # @param d_y: Known length of the y direction aperture width. If None, the function will solve for it.
    # @return: The computed value of the missing variable(s)
    # @throws ValueError: If insufficient information is provided to solve for the missing variable.   
    @staticmethod 
    def function21_39(delta_theta_d_x=None, delta_theta_d_y=None, lambda_=None, d_x=None, d_y=None):
        
        # If both coordinates are missing, calculate them
        if delta_theta_d_x is None and delta_theta_d_y is None:
            if lambda_ is not None and d_x is not None and d_y is not None:
                delta_theta_d_x = lambda_ / d_x
                delta_theta_d_y = lambda_ / d_y
                return delta_theta_d_x, delta_theta_d_y
            else:
                raise ValueError("You must provide lambda and d_x to solve for delta_theta_d_x and lambda and d_y to solve for delta_theta_d_y")
            
        # If delta_theta_d_x is missing and delta_theta_d_y is known, calculate delta_theta_d_x
        elif delta_theta_d_x is None and delta_theta_d_y is not None and d_x is not None:
            if lambda_ is not None:
                delta_theta_d_x = lambda_ / d_x
                return delta_theta_d_x, delta_theta_d_y
            else:
                raise ValueError("You must provide lambda and d_x to solve for delta_theta_d_x")
        
        # If delta_theta_d_y is missing and delta_theta_d_x is known, calculate delta_theta_d_y    
        elif delta_theta_d_y is None and delta_theta_d_x is not None and lambda_ is not None and d_y is not None:
            if d_x is None:
                delta_theta_d_y = lambda_ / d_y
                return delta_theta_d_x, delta_theta_d_y
            else:
                raise ValueError("You must provide lambda and d_y to solve for delta_theta_d_y")
         
        # If lambda_ is missing and delta_theta_d_x is known, calculate lambda_    
        elif lambda_ is None and delta_theta_d_x is not None and d_y is None and delta_theta_d_y is None:
            if d_x is not None:
                lambda_ = delta_theta_d_x * d_x
                return lambda_
            else:
                raise ValueError("You must provide d_x to solve for lamda")
        
        # If lambda_ is missing and delta_theta_d_y is known, calculate lambda_    
        elif lambda_ is None and delta_theta_d_y is not None:
            if d_y is not None:
                lambda_ = delta_theta_d_y * d_y
                return lambda_
            else:
                raise ValueError("You must provide d_y to solve for lamda")
            
        # If d_x is missing, calculate d_x    
        elif d_x is None and delta_theta_d_y is None:
            if lambda_ is not None and delta_theta_d_x is not None:
                d_x = lambda_ / delta_theta_d_x
                return d_x
            else:
                raise ValueError("You must provide lambda and delta_theta_d_x to solve for d_x")
            
        # If d_x is missing, calculate d_x
        elif d_y is None:
            if lambda_ is not None and delta_theta_d_y is not None:
                d_y = lambda_ / delta_theta_d_y
                return d_y
            else:
                raise ValueError("You must provide lambda and delta_theta_d_x to solve for d_y")
        else:
            raise ValueError("Invalid combination of missing variables")
    
    # Solve for the missing variable or coordinates in the equations delta_r_d = (1.22 * lambda_ * z) / d and delta_theta_d = (1.22 * lambda_) / d
    # Where delta_r_d, delta_theta_d, lambda, z, and d may be known and one is missing. delta_r_d and delta_theta_d may both be missing if the rest
    # of the variables are known. This function solves for each variable based on the provided known values. If any required variable(s) is missing,
    # the funciton will compute its value using the given equation.
    # @param delta_r_d: the spacial resolution of the given rectangular aperture. If None, the function will solve for it.
    # @param delta_theta_d: the angular resolution of the given rectangular aperture. If None, the function will solve for it.
    # @param lambda_: Known wavelength. If None, the function will solve for it.
    # @param z: Known separation distance between object plane and entrance pupil plane of the imaging system. If None, the function will solve for it.
    # @param d: Known spacial direction. If None, the function will solve for it.
    # @return: The computed value of the missing variable(s)
    # @throws ValueError: If insufficient information is provided to solve for the missing variable.      
    @staticmethod    
    def function21_40(delta_r_d=None, delta_theta_d=None, lambda_=None, z=None, d=None):
        
        # If both coordinates are missing, calculate them
        if delta_r_d is None and delta_theta_d is None:
            if lambda_ is not None and z is not None and d is not None:
                delta_r_d = (1.22 * lambda_ * z) / d
                delta_theta_d = (1.22 * lambda_) / d
                return delta_r_d, delta_theta_d
            else:
                raise ValueError("You must provide lambda, z and d to solve for delta_r_d and lambda and d to solve for delta_theta_d")
        
        # If delta_r_d is missing and delta_theta_d is known, calculate delta_r_d    
        elif delta_r_d is None and delta_theta_d is not None and z is not None:
            if lambda_ is not None and d is not None:
                delta_r_d = (1.22 * lambda_ * z) / d
                return delta_r_d, delta_theta_d
            else:
                raise ValueError("You must provide lambda, z and d to solve for delta_r_d")
        
        # If delta_theta_d is missing and delta_r_d is known, calculate delta_theta_d    
        elif delta_theta_d is None and delta_r_d is not None and z is None:
            if lambda_ is not None and d is not None:
                delta_theta_d = (1.22 * lambda_) / d
                return delta_r_d, delta_theta_d
            else:
                raise ValueError("You must provide lambda and d to solve for delta_theta_d.")
            
        # If d is missing and delta_r_d is known, calculate d    
        elif d is None and delta_r_d is not None:
            if lambda_ is not None and z is not None:
                d = (1.22 * lambda_ * z) / delta_r_d
                return d
            else:
                raise ValueError("You must provide lambda, z and delta_r_d to solve for d")
            
        # If d is missing and delta_theta_d is known, calculate d    
        elif d is None and delta_theta_d is not None:
            if lambda_ is not None:
                d = (1.22 * lambda_) / delta_theta_d
                return d
            else:
                raise ValueError("You must provide lambda and delta_theta_d to solve for d")
        
        # If lambda_ is missing and delta_r_d is known, calculate lambda_    
        elif lambda_ is None and delta_r_d is not None:
            if d is not None and z is not None:
                lambda_ = (delta_r_d * d) / (1.22 * z)
                return lambda_
            else:
                raise ValueError("You must provide d and z and delta_r_d to solve for lambda")
            
        # If lambda_ is missing and delta_theta_d is known, calculate lambda_    
        elif lambda_ is None and delta_theta_d is not None:
            if d is not None:
                lambda_ = (delta_theta_d * d) / 1.22
                return lambda_
            else:
                raise ValueError("You must provide d and delta_theta_d to solve for lambda")
        
        # If z is missing, calculate z    
        elif z is None and delta_r_d is not None and delta_theta_d is not None:
            if lambda_ is not None and d is not None:
                z = (delta_r_d * d) / (1.22 * lambda_)
                return z
            else:
                raise ValueError("You must provide delta_r_d, lamdba_, and d to solve for z")
            
        else:
            raise ValueError("Invalid combination of missing variables")
            
        
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
    