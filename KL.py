# import numpy as np
# from scipy.optimize import minimize
#
#
# class MyClass:
#     def __init__(self, l, alpha):
#         self.l = l
#         self.alpha = alpha
#         self.P_star = np.array([alpha] + [(1 - alpha) / (l - 1)] * (l - 1))
#
#
#     def equality_constraint(self, p):
#         # Ensure the sum of probabilities equals 1
#         return np.sum(p) - 1
#
#     def inequality_constraints(self, p):
#         # Ensure each probability is less than alpha
#         return self.alpha - p
#
#     def objective_function(self, p):
#         """Calculates the Kullback-Leibler divergence between Q and P*."""
#         Q = np.ones(self.l) / self.l  # Uniform distribution Q
#         # Calculate Kullback-Leibler divergence from P* to Q
#         KL = np.sum(self.P_star * np.log(self.P_star / Q))
#         return KL
#
# if __name__ == '__main__':
#     # Example using scipy.optimize.minimize with SLSQP solver
#     num_variables = 10
#     alpha = 0.5
#     my_instance = MyClass(num_variables, alpha)
#     solution = minimize(my_instance.objective_function,  # Objective function
#                         np.ones(num_variables) / num_variables,  # Initial guess for p
#                         method='SLSQP',  # Sequential Least Squares Programming solver
#                         constraints=[{'type': 'eq', 'fun': my_instance.equality_constraint},
#                                      {'type': 'ineq', 'fun': my_instance.inequality_constraints}])
#     optimal_p = solution.x
#     # Verify that the sum of probabilities is 1
#     sum_p = np.sum(optimal_p)
#     print("Sum of probabilities:", sum_p)
#     # Verify that each probability is less than or equal to alpha
#     less_than_alpha = all(p <= alpha for p in optimal_p)
#     print("Each probability is less than or equal to alpha:", less_than_alpha)
#     # Calculate and print the Kullback-Leibler divergence
#     kl_divergence = my_instance.objective_function(optimal_p)
#     print("Kullback-Leibler divergence:", kl_divergence)

# import numpy as np
# from scipy.optimize import minimize
#
#
# class MyClass:
#     def __init__(self, l, alpha):
#         self.l = l
#         self.alpha = alpha
#         self.P_star = np.array([alpha] + [(1 - alpha) / (l - 1)] * (l - 1))
#         self.Q = np.ones(l) / l  # Uniform distribution Q
#
#     def equality_constraint(self, p):
#         return np.sum(p) - 1
#
#     def inequality_constraints(self, p):
#         return self.alpha - p
#
#     def kl_divergence(self, p, P):
#         """Calculates the Kullback-Leibler divergence from P to Q."""
#         return np.sum(P * np.log(P / self.Q))
#
#     def objective_function(self, p):
#         """Objective function for optimization."""
#         return self.kl_divergence(p, p)
#
# if __name__ == '__main__':
#
#     # Example usage
#     num_variables = 10
#     alpha = 0.5
#     my_instance = MyClass(num_variables, alpha)
#     solution = minimize(my_instance.objective_function,
#                         np.ones(num_variables) / num_variables,
#                         method='SLSQP',
#                         constraints=[{'type': 'eq', 'fun': my_instance.equality_constraint},
#                                      {'type': 'ineq', 'fun': my_instance.inequality_constraints}])
#
#     # Calculate D(P || Q) for the optimized P
#     optimal_p = solution.x
#     dpq = my_instance.kl_divergence(optimal_p, optimal_p)
#
#     # Calculate D(P* || Q)
#     dp_star_q = my_instance.kl_divergence(my_instance.P_star, my_instance.P_star)
#
#     print("D(P || Q) for optimized P:", dpq)
#     print("D(P* || Q):", dp_star_q)

import numpy as np
from scipy.optimize import minimize


class MyClass:
    def __init__(self, l, alpha):
        self.l = l
        self.alpha = alpha
        self.P_star = np.array([alpha] + [(1 - alpha) / (l - 1)] * (l - 1))
        self.Q = np.ones(l) / l  # Uniform distribution Q

    def equality_constraint(self, p):
        return np.sum(p) - 1

    def inequality_constraints(self, p):
        return p - self.alpha

    # additional constraint that enforces that at least one probability
    # in the distribution p is not equal to the value 1/l.
    def non_uniform_constraint(self, p):
        # Enforce that at least one probability is not equal to 1/l
        return 1 - np.isclose(p, 1.0 / self.l).astype(int).sum()

    def kl_divergence(self, P):
        return np.sum(P * np.log(P / self.Q))

    def objective_function(self, p):
        return self.kl_divergence(p)

if __name__ == '__main__':

    # Example usage
    l = 10
    alpha = 0.5
    print(f'l={l}, alpha={alpha}')
    my_instance = MyClass(l, alpha)
    initial_guess = np.ones(l) / l
    initial_guess[0] += 0.01  # Slightly modify the initial guess to not be uniform
    initial_guess[-1] -= 0.01
    solution = minimize(my_instance.objective_function,
                        initial_guess,
                        method='SLSQP',
                        constraints=[{'type': 'eq', 'fun': my_instance.equality_constraint},
                                     {'type': 'ineq', 'fun': my_instance.inequality_constraints},
                                     {'type': 'ineq', 'fun': my_instance.non_uniform_constraint}])

    # Calculate D(P || Q) for the optimized P
    optimal_p = solution.x
    dpq = my_instance.kl_divergence(optimal_p)
    print(f'optimal P={optimal_p}')
    # Calculate D(P* || Q)
    print(f'P_star={my_instance.P_star}')
    dp_star_q = my_instance.kl_divergence(my_instance.P_star)

    print("D(P || Q) for optimized P:", dpq)
    print("D(P* || Q):", dp_star_q)
    print(f'D(P || Q) > D(P* || Q)={dpq>dp_star_q}')
