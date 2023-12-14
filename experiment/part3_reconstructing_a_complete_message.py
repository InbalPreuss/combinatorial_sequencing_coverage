from math import ceil

from matplotlib import pyplot as plt
from scipy.stats import norm, binom
import numpy as np
from part2_reconstructing_a_complete_combinatorial_sequence import run_part2_reconstructing_combinatorial_sequence
from part1_reconstructing_a_single_combinatorial_position import run_part1_reconstructing_single_combinatorial_position


class MessageDecoder:
    def __init__(self, n, t, m, b, l, a, delta, eps):
        """
        Initialize the parameters for message reconstruction.

        :param n: Total number of unique building blocks in each position.
        :param t: Required threshold on the number of observed occurrences.
        :param m: Sequence length for each barcode.
        :param b: Number of letters required to be successfully decoded in each barcode.
        :param l: Number of barcodes in the message.
        :param a: Number of barcodes required to be successfully decoded.
        :param delta: Acceptable error threshold.
        """
        self.n = n
        self.t = t
        self.m = m
        self.b = b
        self.l = l
        self.a = a
        self.delta = delta
        self.eps = eps

    def visualize_values_and_prob(self, values, x_label, P_values, y_label, title):
        plt.plot(values, P_values, marker='o')
        plt.axhline(y=np.sqrt(1 - delta), color='r', linestyle='--')
        plt.xlabel(x_label)
        plt.ylabel(y_label)
        plt.title(f'Finding {title}')
        plt.show()

    def find_T(self):
        """
        Find T using iterative search or solving the normal approximation equation.
        """
        T_values = []
        P_X_T_values = []

        T = 1  # Starting with an initial value for T (pi_R)
        while True:
            # Call Part 1 with T to get π(T)
            _, pi_T, _ = run_part1_reconstructing_single_combinatorial_position(self.n, self.t, self.eps, T)
            # TODO: use binom.sf instead fo sum(pmf)
            # Use π(T) as pi_R in Part 2 to get P_succ_single
            P_succ_single = run_part2_reconstructing_combinatorial_sequence(self.m, self.b, pi_T, "normal")

            # TODO: use binom.sf instead fo sum(pmf)
            # Calculate P(X_T >= a)
            P_X_T = sum(binom.pmf(k, self.l, P_succ_single) for k in range(self.a, self.l + 1))

            T_values.append(T)
            P_X_T_values.append(P_X_T)
            print(f'T={T}, part1 pi_T={pi_T}, part2 P_succ_single={P_succ_single}, part3 P_X_T={P_X_T}, np.sqrt(1 - self.delta)={np.sqrt(1 - self.delta)}, P_X_T >= np.sqrt(1 - self.delta)={P_X_T >= np.sqrt(1 - self.delta)}')

            # Check if P(X_T >= a) meets the condition
            if P_X_T >= np.sqrt(1 - self.delta):
                self.visualize_values_and_prob(values=T_values,x_label='T', P_values=P_X_T_values, y_label='P(X_T >= a)', title='T')
                return T
            T += 5  # Increment T

    def calculate_R_all(self, T):
        """
        Calculate R_all using an iterative search.
        """
        R_all_values = []
        P_E_values = []

        R_all = T * l
        while True:
            # Calculate the probability P(E)
            P_E = min(1, ((R_all + 1) ** self.l) * 2 ** (-R_all * self.epsilon(T, R_all)))

            R_all_values.append(R_all)
            P_E_values.append(1 - P_E)

            print(f'R_all={R_all}, P_E={P_E}, 1-P_E={(1-P_E)}, np.sqrt(1 - self.delta)={np.sqrt(1 - self.delta)},1 - P_E >= np.sqrt(1 - self.delta)={1 - P_E >= np.sqrt(1 - self.delta)}')
            if 1 - P_E >= np.sqrt(1 - self.delta):
                self.visualize_values_and_prob(values=R_all_values,x_label='R_all', P_values=P_E_values, y_label='1 - P(E)', title='R_all')
                return R_all
            R_all = ceil(R_all*1.2)  # Increment R_all

    def epsilon(self, T, R_all):
        """
        Calculate the Kullback-Leibler divergence ε(T).
        """
        Q = [1 / self.l] * l
        P = [(T - 1) / R_all] + [(R_all - T + 1) / ((self.l - 1) * R_all)] * (l - 1)
        return sum(p * np.log(p / q) for p, q in zip(P, Q))#TODO: find KL function

    def run_decoder(self):
        T = self.find_T()
        R_all = self.calculate_R_all(T)
        print(f"Calculated T: {T}, R_all: {R_all}")


if __name__ == "__main__":
    # Example input parameters
    n = 5  # Total number of unique building blocks in each position
    t = 3  # Required threshold on the number of observed occurrences
    m = 120  # Sequence length (for each barcode)
    b = 100  # Number of letters required to be successfully decoded in each barcode
    l = 10  # Number of barcodes in the message
    a = ceil(l*0.8)  # Number of barcodes required to be successfully decoded
    eps = 0.01
    delta = 0.1  # Acceptable error threshold

    decoder = MessageDecoder(n, t, m, b, l, a, delta, eps)
    decoder.run_decoder()
