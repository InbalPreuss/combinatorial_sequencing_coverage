from math import ceil
from typing import Tuple, Union, Any, Iterable

from numpy import ndarray
from scipy.special import kl_div
from matplotlib import pyplot as plt
from scipy.stats import norm, binom
import numpy as np
from calculation.part2_reconstructing_a_complete_combinatorial_sequence import ReconstructingCombinatorialSequence


class MessageDecoder:
    def __init__(self, n, t, eps, R, m, b, method, l, a, delta):
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
        self.R=R
        self.method=method

    def visualize_values_and_prob(self, values, x_label, P_values, y_label, title):
        plt.plot(values, P_values, marker='o')
        plt.axhline(y=np.sqrt(1 - self.delta), color='r', linestyle='--')
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
            P_X_T, P_succ_single, pi_T = self.calc_P_X_T(T=T)

            T_values.append(T)
            P_X_T_values.append(P_X_T)
            print(f'T={T}, part1 pi_T={pi_T}, part2 P_succ_single={P_succ_single}, part3 P_X_T={P_X_T}, np.sqrt(1 - self.delta)={np.sqrt(1 - self.delta)}, P_X_T >= np.sqrt(1 - self.delta)={P_X_T >= np.sqrt(1 - self.delta)}')

            # Check if P(X_T >= a) meets the condition
            if P_X_T >= np.sqrt(1 - self.delta):
                self.visualize_values_and_prob(values=T_values,x_label='T', P_values=P_X_T_values, y_label='P(X_T >= a)', title='T')
                return T
            T += 5  # Increment T

    def calc_R_all(self, T):
        """
        Calculate R_all using an iterative search.
        """
        R_all_values = []
        P_E_values = []

        R_all = T * self.l
        while True:
            # Calculate the probability P(E)
            P_E = self.calc_P_E(T=T,R_all=R_all)

            R_all_values.append(R_all)
            P_E_values.append(1 - P_E)

            print(f'R_all={R_all}, P_E={P_E}, 1-P_E={(1-P_E)}, np.sqrt(1 - self.delta)={np.sqrt(1 - self.delta)},1 - P_E >= np.sqrt(1 - self.delta)={1 - P_E >= np.sqrt(1 - self.delta)}, T={T}')
            if 1 - P_E >= np.sqrt(1 - self.delta):
                # self.visualize_values_and_prob(values=R_all_values,x_label='R_all', P_values=P_E_values, y_label='1 - P(E)', title='R_all')
                return R_all
            R_all = ceil(R_all*1.2)  # Increment R_all

    def epsilon(self, T, R_all):
        """
        Calculate the Kullback-Leibler divergence ε(T).
        """
        Q = [1 / self.l] * self.l
        P = [(T - 1) / R_all] + [(R_all - T + 1) / ((self.l - 1) * R_all)] * (self.l - 1)
        return np.sum(kl_div(P, Q))

    def run_decoder(self):
        T = self.find_T()
        R_all = self.calc_R_all(T)
        print(f"Calculated T: {T}, R_all: {R_all}")

    def calc_P_X_T(self, T: int) -> tuple[Any, float, Any]:
        # Use π(T) as pi_R in Part 2 to get P_succ_single
        reconstructor = ReconstructingCombinatorialSequence(n=self.n, t=self.t, eps=self.eps, R=T, m=self.m, b=self.b, method=self.method)
        P_succ_single, pi_T = reconstructor.calculate_probability()

        # Calculate P(X_T >= a)
        P_X_T = binom.sf(self.a - 1, self.l, P_succ_single)

        return P_X_T, P_succ_single, pi_T

    def calc_P_E(self, T, R_all):
        try:
            result = min(1, ((R_all + 1) ** self.l) * 2 ** (-R_all * self.epsilon(T, R_all)))
        except OverflowError:
            # Handle the overflow, for example, set result to a default value
            result = 1

        return result


if __name__ == "__main__":
    ##########
    # Part 1 #
    ##########
    n = 5  # Total number of unique building blocks in each position
    t = 3  # Required threshold on the number of observed occurrences
    eps = 0.01
    R = 40  # Acceptable error threshold

    ##########
    # Part 2 #
    ##########
    m = 120  # Total number of letters in the sequence
    b = 100  # Number of letters required to be successfully reconstructed
    method = "binomial"

    ##########
    # Part 3 #
    ##########
    # Example input parameters
    # n = 5  # Total number of unique building blocks in each position
    # t = 3  # Required threshold on the number of observed occurrences
    # m = 120  # Sequence length (for each barcode)
    # b = 100  # Number of letters required to be successfully decoded in each barcode
    l = 10  # Number of barcodes in the message
    a = ceil(l*0.8)  # Number of barcodes required to be successfully decoded
    # eps = 0.01
    delta = 0.1  # Acceptable error threshold

    decoder = MessageDecoder(n=n, t=t, eps=eps, R=R, m=m, b=b, method=method, l=l, a=a, delta=delta)
    decoder.run_decoder()
