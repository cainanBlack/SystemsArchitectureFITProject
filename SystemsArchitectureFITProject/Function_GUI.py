import tkinter as tk
from tkinter import messagebox
from SystemsArchitectureFITProject.Functions.Functions import Functions


# Assuming the Functions class with your existing functions is imported

class FunctionSolverApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Function Solver")

        # Create the UI elements
        self.create_widgets()

    def create_widgets(self):
        # Function to solve for missing variables in 21_18
        self.func21_18_button = tk.Button(self.root, text="Solve Function 21_18", command=self.solve_function21_18)
        self.func21_18_button.grid(row=0, column=0, padx=250, pady=10)

        # Function to solve for missing variables in 21_30
        self.func21_30_button = tk.Button(self.root, text="Solve Function 21_30", command=self.solve_function21_30)
        self.func21_30_button.grid(row=1, column=0, padx=250, pady=10)

        # Function to solve for missing variables in 21_44
        self.func21_44_button = tk.Button(self.root, text="Solve Function 21_44", command=self.solve_function21_44)
        self.func21_44_button.grid(row=2, column=0, padx=250, pady=10)

        # Function to solve for missing variables in 21_58
        self.func21_58_button = tk.Button(self.root, text="Solve Function 21_58", command=self.solve_function21_58)
        self.func21_58_button.grid(row=3, column=0, padx=250, pady=10)

        # Function to solve for missing variables in 21_59
        self.func21_59_button = tk.Button(self.root, text="Solve Function 21_59", command=self.solve_function21_59)
        self.func21_59_button.grid(row=4, column=0, padx=250, pady=10)

        # Function to solve for missing variables in 21_85
        self.func21_85_button = tk.Button(self.root, text="Solve Function 21_85", command=self.solve_function21_85)
        self.func21_85_button.grid(row=5, column=0, padx=250, pady=10)

        # Function to solve for missing variables in 21_104
        self.func21_104_button = tk.Button(self.root, text="Solve Function 21_104", command=self.solve_function21_104)
        self.func21_104_button.grid(row=6, column=0, padx=250, pady=10)

        # Function to solve for missing variables in 21_63
        self.func21_63_button = tk.Button(self.root, text="Solve Function 21_63", command=self.solve_function21_63)
        self.func21_63_button.grid(row=7, column=0, padx=250, pady=10)

    def solve_function21_18(self):
        # Create input dialog for function 21_18
        self.open_input_dialog("function21_18")

    def solve_function21_30(self):
        # Create input dialog for function 21_30
        self.open_input_dialog("function21_30")

    def solve_function21_44(self):
        # Create input dialog for function 21_44
        self.open_input_dialog("function21_44")

    def solve_function21_58(self):
        # Create input dialog for function 21_58
        self.open_input_dialog("function21_58")

    def solve_function21_59(self):
        # Create input dialog for function 21_59
        self.open_input_dialog("function21_59")

    def solve_function21_85(self):
        # Create input dialog for function 21_85
        self.open_input_dialog("function21_85")

    def solve_function21_104(self):
        # Create input dialog for function 21_104
        self.open_input_dialog("function21_104")

    def solve_function21_63(self):
        # Create input dialog for function 21_63
        self.open_input_dialog("function21_63")

    def open_input_dialog(self, func_name):
        # Open dialog box to input parameters for the given function
        def get_values():
            values = {key: entry.get() for key, entry in entry_widgets.items()}
        
            # Explicitly replace 'lambda' with 'lambda_param'
            values = {('lambda_param' if key == 'lambda' else key): (float(val) if val else None) for key, val in values.items()}
        
            try:
                # Call the function with the corrected argument names
                if func_name == "function21_18":
                    result = Functions.function21_18(**values)
                elif func_name == "function21_30":
                    result = Functions.function21_30(**values)
                elif func_name == "function21_44":
                    result = Functions.function21_44(**values)
                elif func_name == "function21_58":
                    result = Functions.function21_58(**values)
                elif func_name == "function21_59":
                    result = Functions.function21_59(**values)
                elif func_name == "function21_85":
                    result = Functions.function21_85(**values)
                elif func_name == "function21_104":
                    result = Functions.function21_104(**values)
                elif func_name == "function21_63":
                    result = Functions.function_21_63(**values)
        
                # Show the result
                messagebox.showinfo("Result", f"The result is: {result}")
            except ValueError as e:
                messagebox.showerror("Error", str(e))

        # Create a window for the input
        input_window = tk.Toplevel(self.root)
        input_window.title(f"Input for {func_name}")

        # Create a dictionary of entry widgets for the function
        entry_widgets = {}

        # Generate the inputs based on the function
        if func_name == "function21_18":
            labels = ["f", "x", "lambda", "d_i"]
        elif func_name == "function21_30":
            labels = ["I_f", "O_f", "H_f"]
        elif func_name == "function21_44":
            labels = ["n_1", "P", "T"]
        elif func_name == "function21_58":
            labels = ["k", "r", "gamma_n"]
        elif func_name == "function21_59":
            labels = ["phi", "gamma", "k"]
        elif func_name == "function21_85":
            labels = ["phi", "a", "rho", "theta", "Z"]
        elif func_name == "function21_104":
            labels = ["h_f", "o_n_f", "k_bar", "var_h_f", "p_sigma_n_squared"]
        elif func_name == "function21_63":
            labels = ["Gamma_p", "Delta_x", "r0"]

        for i, label in enumerate(labels):
            tk.Label(input_window, text=label).grid(row=i, column=0, padx=150, pady=5)
            entry = tk.Entry(input_window)
            entry.grid(row=i, column=1, padx=150, pady=5)
            entry_widgets[label] = entry

        # Add a button to submit the input and call the respective function
        submit_button = tk.Button(input_window, text="Submit", command=get_values)
        submit_button.grid(row=len(labels), columnspan=2, pady=10)

if __name__ == "__main__":
    root = tk.Tk()
    app = FunctionSolverApp(root)
    root.mainloop()
