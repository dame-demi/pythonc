# Import the function from your module
from KMEquation import solve


def main():
    # Example input values
    initial_conditions = (0.0215e-6, 0.1)  # Example initial conditions
    rb_initial = 0.0215e-6
    CL = 1300.0
    Pv = 1160000
    p_ambient1 = 116000
    p_ambient2 = 116000
    sigma = 0.25
    n = 1.0
    mu = 0.24e-3
    kappa = 8.0e-6
    rho = 595.59
    py_list = [0.0, 1.5e-9]  # Example list
    specific_time = 1.5e-9

    # Call the solve function
    R_values = solve(
        initial_conditions, rb_initial, CL, Pv, p_ambient1, p_ambient2,
        sigma, n, mu, kappa, rho, py_list, specific_time
    )

    # Print the R_values list and its type to ensure it is correct
    print(f"R_values: {R_values}, type: {type(R_values)}")


