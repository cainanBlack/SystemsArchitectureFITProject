
from SystemsArchitectureFITProject.Functions.Functions_for_Functions import FunctionsFor21_19, FunctionsFor21_30, FunctionsFor21_32, FunctionsFor21_58, FunctionsFor21_59, FunctionsFor21_85
import math
from scipy.optimize import fsolve, root
import numpy as np
from mpmath import mp
from scipy.integrate import quad
from scipy.optimize import root_scalar
from scipy.special import jn  # Correctly importing the Bessel function

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
    def function21_18(f=None, x=None, lambda_param=None, d_i=None):
        # If f is missing, calculate f (average frequency)
        if f is None:
            if x is not None and lambda_param is not None and d_i is not None:
                f = x / (lambda_param * d_i)
                return f
            else:
                raise ValueError("You must provide x, lambda, and d_i to solve for f.")
        
        # If x is missing, calculate x (average value of x)
        elif x is None:
            if f is not None and lambda_param is not None and d_i is not None:
                x = f * lambda_param * d_i
                return x
            else:
                raise ValueError("You must provide f, lambda, and d_i to solve for x.")
        
        # If lambda is missing, calculate lambda
        elif lambda_param is None:
            if f is not None and x is not None and d_i is not None:
                lambda_param = x / (f * d_i)
                return lambda_param
            else:
                raise ValueError("You must provide f, x, and d_i to solve for lambda.")
        
        # If d_i is missing, calculate d_i
        elif d_i is None:
            if f is not None and x is not None and lambda_param is not None:
                d_i = x / (f * lambda_param)
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
    
    # Solve for any missing variable in the equation:
    # h(r) = (D / (2 * lamda * d_i))^2 * J1((2 * pi * r * D) / (2 * lamda * d_i))
    # Arguments:
    # @param h     The value of h(r), can be None if solving for h.
    # @param r     The value of radius r, can be None if solving for r.
    # @param D     The value of diameter D, can be None if solving for D.
    # @param lamda The wavelength (lambda), can be None if solving for lambda.
    # @param d_i The diameter of the disk, can be None if solving for d_i.
    # Returns:
    # The solved value for the unknown variable.
    @staticmethod
    def function21_29(h=None, r=None, D=None, lamda=None, d_i=None):
        # Check if any of the required variables are missing
        if (h is None and r is None) or (h is None and D is None) or (h is None and lamda is None) or (h is None and d_i is None):
            raise ValueError("At least one variable must be provided")
    
        # Define the equation in terms of the unknown variable
        def equation_to_solve(var_value):
            if h is not None and r is None:
                # Solve for r
                equation = (D / (2 * lamda * d_i))**2 * jn(1, (2 * np.pi * var_value * D) / (2 * lamda * d_i)) - h
            elif h is not None and D is None:
                # Solve for D
                equation = (var_value / (2 * lamda * d_i))**2 * jn(1, (2 * np.pi * r * var_value) / (2 * lamda * d_i)) - h
            elif h is not None and lamda is None:
                # Solve for lamda
                equation = (D / (2 * var_value * d_i))**2 * jn(1, (2 * np.pi * r * D) / (2 * var_value * d_i)) - h
            elif h is not None and d_i is None:
                # Solve for d_i
                equation = (D / (2 * lamda * var_value))**2 * jn(1, (2 * np.pi * r * D) / (2 * lamda * var_value)) - h
            elif r is not None and h is None:
                # Solve for h
                equation = (D / (2 * lamda * d_i))**2 * jn(1, (2 * np.pi * r * D) / (2 * lamda * d_i)) - var_value
            return equation
    
        # Initial guess for the unknown variable (use current known value for initial guess)
        if r is not None:
            initial_guess = r
        elif D is not None:
            initial_guess = D
        elif lamda is not None:
            initial_guess = lamda
        elif d_i is not None:
            initial_guess = d_i
        else:
            raise ValueError("At least one variable must be provided")
    
        # Solve the equation using fsolve
        solved_value = fsolve(equation_to_solve, initial_guess)
        return solved_value[0]

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

    # Solve for the optical transfer function equation 
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
    def function21_32(Hf = None, fx=None, fy=None, lambda_=None, di=None, Dx=None, Dy=None):
        # If  Hf is missing, calculate  Hf
        if Hf is None:
            if fx is not None and fy is not None and lambda_ is not None and di is not None and  Dx is not None and  Dy is not None:

                tri_x = FunctionsFor21_32.triangular((fx * lambda_ * di) / Dx)
                tri_y = FunctionsFor21_32.triangular((fy * lambda_ * di) / Dy)
                
                return  tri_x * tri_y
            
            # If lamba is missing, calculate lamba
            elif lambda_ is None and fx is not None and fy is not None and di is not None and Dx is not None and Dy is not None:
                return FunctionsFor21_32.compute_Lambda(fx, fy, di, Dx, Dy)
            
            else:
                raise ValueError("You must provide fx, fy, lambda, di, Dx, Dy to solve for H_f.")

        elif Hf is not None:

            # if Dx is None and Dy is None:
            #     raise ValueError("You must provide Dx or Dy to solve any missing pieace.")

            # If di is missing, calculate Dx
            if di is None and fx is not None and fy is not None and lambda_ is not None and Dx is not None and Dy is not None :
                return FunctionsFor21_32.compute_di(Hf, fx, fy, lambda_, Dx, Dy)

            elif di is not None:
                # If Dx is missing, calculate Dx
                if Dx is None and fx is not None and lambda_ is not None and fy is None and Dy is None: 
                    return FunctionsFor21_32.compute_Dx(Hf, fx, lambda_, di)
                
                # If Dy is missing, calculate Dy
                elif Dy is None and fy is not None and lambda_ is not None and Dx is None and fx is None :
                    return FunctionsFor21_32.compute_Dy(Hf, fy, lambda_, di)

                # If fx is missing, calculate fx
                elif fx is None and Dx is not None and lambda_ is not None and Dy is None and fy is None :
                    return FunctionsFor21_32.compute_fx(Hf, fy, lambda_, Dx)
                
                # If fy is missing, calculate fy
                elif fy is None and Dy is not None and lambda_ is not None and (Dx is None or fx is None) :
                    return FunctionsFor21_32.compute_fy(Hf, fy, lambda_, Dx)

                #at this point what is left to be solved                
                else:
                    raise ValueError("You must provide di to solve any missing pieace.")

            else:
                raise ValueError("You must provide di to solve any missing pieace.")


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

    # Solves the Kolmogorov Power Spectral Density (PSD) equation for a missing variable or calculates PSD directly
    # @param missing (str): The variable to solve for ('Phi', 'C_n', 'K'), or None if no variable is missing
    # @param Phi (float): The power spectral density Φ_n^K (optional if solving for Phi)
    # @param C_n (float): The refractive index structure constant (optional if solving for C_n)
    # @param K (float): The wavenumber (optional if solving for K)
    # @return float: The value of the missing variable
    @staticmethod
    def function21_45(Phi=None, C_n=None, K=None):

        # Solve for Φ_n^K (Phi)
        if Phi is None:
            if C_n is not None or K is not None:
                return 0.033 * (C_n**2) * (K**(-11/3))
            else:
                raise ValueError("C_n and K must be provided to solve for Phi.")
        
        # Solve for C_n
        elif C_n is None:
            if Phi is not None or K is not None:
                return (Phi / (0.033 * K**(-11/3)))**0.5
            else:
                raise ValueError("Phi and K must be provided to solve for C_n.")
        
        # Solve for K
        elif K is None:
            if Phi is not None or C_n is not None:
                return (Phi / (0.033 * C_n**2))**(-3/11)
            else:
                raise ValueError("Phi and C_n must be provided to solve for K.")
        
        else:
            raise ValueError("Must be missing a variable to solve for.")

    # Solves the Kolmogorov Power Spectral Density (V-variant) Φ_n^V(K) equation for a missing variable.
    # @param Phi (float): Power spectral density Φ_n^V(K) (optional if solving for Phi)
    # @param C_n (float): Refractive index structure constant (optional if solving for C_n)
    # @param K (float): Wavenumber (Constant not to be solved for)
    # @param K_0 (float): Outer scale of turbulence (Constant not to be solved for)
    # @param K_m (float): Cutoff wavenumber (Constant not to be solved for)
    # @return float: The calculated value of the missing variable or the result if no variable is missing
    @staticmethod
    def function21_48(Phi=None, C_n=None, K=None, K_0=None, K_m=None):
        
        if Phi is None:
            if C_n is not None and K is not None and K_0 is not None and K_m is not None:
                numerator = 0.33 * C_n**2
                denominator = (K**2 + K_0**2)**(11*6)
                exponential = math.exp(-((K**2)/(K_m**2)))
                Phi = ((numerator * exponential) / denominator)
                return Phi
            else:
                raise ValueError("Not Valid 1") 
        
        elif C_n is None:
            if Phi is not None and K is not None and K_0 is not None and K_m is not None:
                C_n = math.sqrt((Phi * (K**2 + K_0**2)**11/6)/ math.exp(-(K**2/K_m**2) * 0.033))
                return C_n
            else:
                raise ValueError("Not Valid 2")
            
        elif K or K_0 or K_m is not None:
            raise NotImplementedError("Constant provided by another equation not solved for here.")           
                
        else:
            raise ValueError("Invalid variable to solve for.")
        
    # Solve for Cn^2(h) or find the missing parameter (h, v, or A) in the given formula.
    # @param h (float): Height in meters. If None, the function will solve for h.
    # @param v (float): Frequency. If None, the function will solve for v.
    # @param A (float): Constant coefficient. If None, the function will solve for A.
    # @param Cn2_value (float): The value of Cn^2(h). Required if solving for a missing parameter.
    # @returns float: The solved or computed parameter.
    @staticmethod
    def function21_52(h=None, v=None, A=None, Cn2_value=None):
        def equation_to_solve(x, missing_param):
            if missing_param == 'h':
                return (
                    5.94e-53 * (v / 27)**2 * x**10 * np.exp(-x / 1000)
                    + 2.7e-16 * np.exp(-x / 1500)
                    + A * np.exp(-x / 100)
                    - Cn2_value
                )
            elif missing_param == 'v':
                return (
                    5.94e-53 * (x / 27)**2 * h**10 * np.exp(-h / 1000)
                    + 2.7e-16 * np.exp(-h / 1500)
                    + A * np.exp(-h / 100)
                    - Cn2_value
                )
            elif missing_param == 'A':
                return (
                    5.94e-53 * (v / 27)**2 * h**10 * np.exp(-h / 1000)
                    + 2.7e-16 * np.exp(-h / 1500)
                    + x * np.exp(-h / 100)
                    - Cn2_value
                )
            else:
                raise ValueError("Invalid missing parameter.")

        # Check which parameter is missing and solve for it
        if h is None:
            if Cn2_value is None or v is None or A is None:
                raise ValueError("Missing inputs: Cn2_value, v, and A must be provided to solve for h.")
            missing_param = 'h'
            initial_guess = 1000  # Initial guess for h
        elif v is None:
            if Cn2_value is None or h is None or A is None:
                raise ValueError("Missing inputs: Cn2_value, h, and A must be provided to solve for v.")
            missing_param = 'v'
            initial_guess = 30  # Initial guess for v
        elif A is None:
            if Cn2_value is None or h is None or v is None:
                raise ValueError("Missing inputs: Cn2_value, h, and v must be provided to solve for A.")
            missing_param = 'A'
            initial_guess = 1e-15  # Initial guess for A
        elif Cn2_value is None:
            if A is None or h is None or v is None:
                raise ValueError("Missing inputs: Cn2_value, h, and v must be provided to solve for A.")
            # If all inputs are provided, calculate Cn^2(h)
            term1 = 5.94e-53 * (v / 27)**2 * h**10 * np.exp(-h / 1000)
            term2 = 2.7e-16 * np.exp(-h / 1500)
            term3 = A * np.exp(-h / 100)
            return term1 + term2 + term3
        else:
            raise ValueError("Must leave one parameter empty to solve")

        # Solve for the missing parameter
        solved_value = fsolve(equation_to_solve, initial_guess, args=(missing_param))
        return solved_value[0]
        
    # Solves for the missing variable: k, r, or gamma_n.
    # Equation Components:
    #     𝛤𝑛(𝑟) = ∫ 𝑑𝑘⃗ 𝛷𝑛(𝑘⃗)𝑒(−𝑗𝑘⃗ ∗𝑟)
    # @param k (tuple) 𝛷𝑛(𝑘⃗)
    # @param r (list) The position vector in real space
    # @param gamma_n (complex) The value of Γₙ(r)
    # @return: The computed value of the equation
    @staticmethod
    def function21_58(r=None, gamma_n=None, k=None):

        def phi_n(k):
            # Example Φₙ(k⃗): Gaussian form
            return np.exp(-np.linalg.norm(k)**2)

        if gamma_n is None:
            if k is not None and r is not None:
                if type(k) is not tuple:
                    raise ValueError("Incorrect type for K! Entered: ", type(k), ", Required: tuple")
                if type(r) is not list:
                    raise ValueError("Incorrect type for R! Entered: ", type(r), " Required: list")
                
                return FunctionsFor21_58.solveForGamma_N(phi_n, r)
            else:
                raise ValueError("phi_n and r are required to solve for gamma_n.")

        elif k is None:
            if r is not None and gamma_n is not None:
                if type(gamma_n) is not complex and type(gamma_n) is not int and type(gamma_n) is not float:
                    raise ValueError("Incorrect type for gamma_n! Entered: ", type(gamma_n), ", Required: complex, int, or float")
                if type(r) is not list:
                    raise ValueError("Incorrect type for R! Entered: ", type(r), " Required: list")
                
                return FunctionsFor21_58.solveForPhi_n(gamma_n, r)
            else:
                raise ValueError("r and gamma_n are required to solve for phi_n.")

        elif r is None:
            if k is not None and gamma_n is not None:
                if type(k) is not tuple:
                    raise ValueError("Incorrect type for K! Entered: ", type(k), ", Required: tuple")
                if type(gamma_n) is not complex and type(gamma_n) is not int and type(gamma_n) is not float:
                    raise ValueError("Incorrect type for gamma_n! Entered: ", type(gamma_n), ", Required: complex, int, or float")
                
                return FunctionsFor21_58.solveForR(phi_n, gamma_n)
            else:
                raise ValueError("phi_n and gamma_n are required to solve for r.")
            
        else:
            raise ValueError("Invalid missing variable. Choose 'gamma_n', 'k', or 'r'.")

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
            if isinstance(k, float):
                k = (k, k, k)  # Make k a tuple if it's a float
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
            if isinstance(k, float):
                k = (k, k, k)  # Make k a tuple if it's a float
            x, y, z = k  # Now k is guaranteed to be unpackable
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
        
    # Solves for the missing variable ('Gamma_p', 'Delta_x', or 'r0') in the equation:
    #     Gamma_p = exp(-6.88 * Delta_x^(5/3) / (2 * r0))
    # @param missing_var The variable to solve for ('Gamma_p', 'Delta_x', or 'r0')
    # @param Gamma_p The given Gamma_p value (optional)
    # @param Delta_: The given Delta_x value (optional)
    # @param r0 The given r0 value (optional)
    # @return The solved value of the missing variable.
    @staticmethod
    def function21_63(Gamma_p=None, Delta_x=None, r0=None):
        mp.dps = 10  # Set precision level for mpmath

        # Check if solving for Gamma_p (if Delta_x and r0 are given)
        if Gamma_p is None:
            if Delta_x is not None and r0 is not None:
                # Direct calculation of Gamma_p (ensure result is real)
                return float(mp.exp(-6.88 * Delta_x**(5/3) / (2 * r0)))  # Convert to float
            
            else:
                raise ValueError("Delta_x and r0 are required to solve for Gamma_p.")
        
        # Check if solving for Delta_x (if Gamma_p and r0 are given)
        elif Delta_x is None:
            if Gamma_p is not None and r0 is not None:
                def equation(vars):
                    # Equation to solve for Delta_x
                    result = (mp.log(Gamma_p) * -2 * r0 / 6.88)**(3/5) - vars[0]
                    return [float(result.real)]  # Ensure result is real and convert to float
                
                # Initial guess for Delta_x (if not provided)
                initial_guess = [Delta_x if Delta_x is not None else 1.0]
                solution = root(equation, initial_guess, method='hybr', options={'xtol': 1e-12, 'maxfev': 100000})
                return solution.x[0]  # Return the first solution
                
            else:
                raise ValueError("Gamma_p and r0 are required to solve for Delta_x.")
        
        # Check if solving for r0 (if Gamma_p and Delta_x are given)
        elif r0 is None:
            if Gamma_p is not None and Delta_x is not None:
                def equation(vars):
                    # Equation to solve for r0
                    result = (mp.log(Gamma_p) * -6.88 * Delta_x**(5/3) / 2) - vars[0]
                    return [float(result.real)]  # Ensure result is real and convert to float
                
                # Initial guess for r0 (if not provided)
                initial_guess = [r0 if r0 is not None else 1.0]
                solution = root(equation, initial_guess, method='hybr', options={'xtol': 1e-12, 'maxfev': 100000})
                return solution.x[0]  # Return the first solution

            else:
                raise ValueError("Gamma_p and Delta_x are required to solve for r0.")
        
        else:
            raise ValueError("No variable is missing. Please leave one of Gamma_p, Delta_x, or r0 as None.")

    # Solves for a missing variable in the Fried parameter equation: r_new = (lambda_new / lambda_old)^(6/5) * r_old
    # @param missing (str): The variable to solve for ('r_new', 'r_old', 'lambda_new', or 'lambda_old')
    # @param r_new (float): The new Fried parameter (optional if solving for r_new)
    # @param r_old (float): The old Fried parameter (optional if solving for r_old)
    # @param lambda_new (float): The new wavelength (optional if solving for lambda_new)
    # @param lambda_old (float): The old wavelength (optional if solving for lambda_old)
    # @return float: The calculated value of the missing variable
    @staticmethod
<<<<<<< Updated upstream
    def function21_68(missing, r_new=None, r_old=None, lambda_new=None, lambda_old=None): #TODO NOTE REMOVE MISSING VALUE PLZ
        
        # Validate input
        if missing not in ['r_new', 'r_old', 'lambda_new', 'lambda_old']:
            raise ValueError("Invalid variable to solve for. Choose 'r_new', 'r_old', 'lambda_new', or 'lambda_old'.")
        
=======
    def function21_68(r_new=None, r_old=None, lambda_new=None, lambda_old=None):
         
>>>>>>> Stashed changes
        # Solve for r_new
        if r_new is None:
            if r_old is not None or lambda_new is not None or lambda_old is not None:
                return (lambda_new / lambda_old)**(6/5) * r_old
            else:
                raise ValueError("r_old, lambda_new, and lambda_old must be provided to solve for r_new.")
        
        # Solve for r_old
        elif r_old is None:
            if r_new is not None or lambda_new is not None or lambda_old is not None:
                return r_new / (lambda_new / lambda_old)**(6/5)
            else:
                raise ValueError("r_new, lambda_new, and lambda_old must be provided to solve for r_old.")
        
        # Solve for lambda_new
        elif lambda_new is None:
            if r_new is not None or r_old is not None or lambda_old is not None:
                return lambda_old * (r_new / r_old)**(5/6)
            else:
                raise ValueError("r_new, r_old, and lambda_old must be provided to solve for lambda_new.")
        
        # Solve for lambda_old
        elif lambda_old is None:
            if r_new is not None or r_old is not None or not lambda_new is None:
                return lambda_new / (r_new / r_old)**(5/6)
            else:
                raise ValueError("r_new, r_old, and lambda_new must be provided to solve for lambda_old.")
            
        else:
            raise ValueError("Invalid combination of missing variables")

    # Solves for the missing value (theta_0, lambda_value, or L) in the given equation.
    # @param theta_0 (float): The known value of θ₀ (must be in range of 0.1e-3 to 0.9e-12).
    # @param lambda_value (float): The wavelength λ (must be in range of 0.1e-3 to 0.9e-8).
    # @param L (float): The upper limit of the integral (must be in range of 1 to 13719818).
    # @param cn_squared_func (function): A function Cn²(z) that describes the refractive index structure parameter as a function of z.
    @staticmethod
    def function21_74(theta_0=None, lambda_value=None, L=None):    
        def cn_squared_func(z):
            return 1e-14  # Example value, adjust as needed 

        # Define the integral
        def calculate_integral(L_val):
            def integrand(z):
                return cn_squared_func(z) * (z ** (5 / 3))
            result, _ = quad(integrand, 0, L_val)
            return result

        # If solving for theta_0
        if theta_0 is None:
            if lambda_value is None or L is None:
                raise ValueError("Both λ and L must be provided to solve for θ₀.")
            integral_result = calculate_integral(L)
            theta_0 = (58.1 * 10**-3 * lambda_value**(6/5)) * (integral_result**(-3/5))
            return theta_0

        # If solving for lambda_value
        if lambda_value is None:
            if theta_0 is None or L is None:
                raise ValueError("Both θ₀ and L must be provided to solve for λ.")
            integral_result = calculate_integral(L)

            def lambda_equation(lam):
                return theta_0 - (58.1 * 10**-3 * lam**(6/5)) * (integral_result**(-3/5))

            result = root_scalar(lambda_equation, bracket=[1e-10, 1], method='bisect')
            return result.root

        # If solving for L
        if L is None:
            if theta_0 is None or lambda_value is None:
                raise ValueError("Both θ₀ and λ must be provided to solve for L.")

            def L_equation(L_val):
                integral_result = calculate_integral(L_val)
                return theta_0 - (58.1 * 10**-3 * lambda_value**(6/5)) * (integral_result**(-3/5))

            result = root_scalar(L_equation, bracket=[1e-2, 1e6], method='bisect')
            return result.root

        # If none of the above, raise an error
        raise ValueError("One of θ₀, λ, or L must be set to None (missing).")

    # Solves for the missing variable in the equation:
    #     I_SE = O * H_SE * H_0
    # This version works with scalar (single number) inputs only.
    # @param I_SE: Output signal (scalar or None if missing).
    # @param O: Input signal (scalar or None if missing).
    # @param H_SE: Transfer function H_SE (scalar or None if missing).
    # @param H_0: Transfer function H_0 (scalar or None if missing).
    # @return: The value of the missing variable as a single number.
    # @throws ValueError: If more than one variable is missing or none is missing.
    @staticmethod
    def function21_78(I_SE=None, O=None, H_SE=None, H_0=None):
        # Check how many variables are missing
        missing_count = sum(1 for v in [I_SE, O, H_SE, H_0] if v is None)
        if missing_count != 1:
            raise ValueError("Only one variable can be missing.")
    
        # Solve for the missing variable
        if I_SE is None:
            return O * H_SE * H_0
        elif O is None:
            return I_SE / (H_SE * H_0)
        elif H_SE is None:
            return I_SE / (O * H_0)
        elif H_0 is None:
            return I_SE / (O * H_SE)      

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
