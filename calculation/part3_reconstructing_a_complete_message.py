import utils
from calculation.part2_reconstructing_a_complete_combinatorial_sequence import ReconstructingCombinatorialSequence
import utils as uts

from math import ceil
from typing import Tuple, Union, Any, Iterable
import scipy as scp
from numpy import ndarray
from scipy.special import kl_div
from matplotlib import pyplot as plt
from scipy.stats import norm, binom, multinomial
import numpy as np


DIR_PLOT_PATH = 'plots/'
DIR_PLOT_PATH_PART_3 = f'{DIR_PLOT_PATH}calculation_part3_plots/'



class MessageDecoder:
    def __init__(self, n, t, eps, R, m, b, method, l, a, delta, P_E_method):
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
        self.R = R
        self.method = method
        self.P_E_method = P_E_method

    def visualize_values_and_prob(self, values, x_label, P_values, y_label, title, is_save_plot=False, is_part_3=False):
        plt.plot(values, P_values, marker='o')
        if is_part_3:
            plt.axhline(y=np.sqrt(1 - self.delta), color='r', linestyle='--')
            plt.text(values[-1]+0.5,np.sqrt(1 - self.delta), r'$\sqrt{1-\delta}$', color='red', verticalalignment='bottom',
                     fontsize=12)

        plt.xlabel(x_label, fontsize=20)
        plt.ylabel(y_label, fontsize=20)
        plt.yticks(fontsize=15)
        plt.xticks(fontsize=15)
        plt.ylim(-0.1,1.1)
        plt.subplots_adjust(bottom=0.2)
        plt.title(f'Finding {title}')
        if is_save_plot:
            uts.make_dir(DIR_PLOT_PATH)
            uts.make_dir(DIR_PLOT_PATH_PART_3)
            plt.savefig(f"{DIR_PLOT_PATH_PART_3}/{title}_prob_with_different_R_recover_single_letter_part1_n={n},t={t},R={R}.svg",
                        format='svg')

        plt.show()

    def find_T(self):
        """
        Find T using iterative search or solving the normal approximation equation.
        """
        T_values = []
        P_X_T_values = []
        P_succ_single_values = []
        pi_T_values = []

        T = 1  # Starting with an initial value for T (pi_R)
        while True:
            P_X_T, P_succ_single, pi_T = self.calc_P_X_T(T=T)

            T_values.append(T)
            P_X_T_values.append(P_X_T)
            P_succ_single_values.append(P_succ_single)
            pi_T_values.append(pi_T)
            print(
                f'T={T}, part1 pi_T={pi_T}, part2 P_succ_single={P_succ_single}, part3 P_X_T={P_X_T}, np.sqrt(1 - self.delta)={np.sqrt(1 - self.delta)}, P_X_T >= np.sqrt(1 - self.delta)={P_X_T >= np.sqrt(1 - self.delta)}')

            # Check if P(X_T >= a) meets the condition
            if P_X_T >= np.sqrt(1 - self.delta):
                self.visualize_values_and_prob(values=T_values, x_label='T', P_values=pi_T_values,
                                               y_label='P', title='Part 1 P', is_save_plot=True)
                self.visualize_values_and_prob(values=T_values, x_label='T', P_values=P_succ_single_values,
                                               y_label='P', title='Part 2 P', is_save_plot=True)
                self.visualize_values_and_prob(values=T_values, x_label='T', P_values=P_X_T_values,
                                               y_label='P(($X_{\mathrm{T}}$) >= a)', title='Part 3 P', is_save_plot=True, is_part_3=True)
                # self.visualize_values_and_prob(values=T_values, x_label='T', P_values=P_X_T_values,
                #                                y_label='P(X_T >= a)', title='T')
                return T
            T += 1  # Increment T

    def calc_R_all(self, T):
        """
        Calculate R_all using an iterative search.
        """
        R_all_values = []
        P_E_values = []

        R_all = T * self.l
        while True:
            # Calculate the probability P(E)
            if P_E_method == 'simulation':
                P_E = self.simu_P_E(T=T, R_all=R_all)
            else:
                P_E = self.calc_P_E(T=T, R_all=R_all)

            R_all_values.append(R_all)
            P_E_values.append(1 - P_E)

            print(
                f'R_all={R_all}, P_E={P_E}, 1-P_E={(1 - P_E)}, np.sqrt(1 - self.delta)={np.sqrt(1 - self.delta)},1 - P_E >= np.sqrt(1 - self.delta)={1 - P_E >= np.sqrt(1 - self.delta)}, T={T}')
            if 1 - P_E >= np.sqrt(1 - self.delta):
                self.visualize_values_and_prob(values=R_all_values,x_label='$R_{\mathrm{all}}$', P_values=P_E_values, y_label='1 - P(E)', title='R_all', is_save_plot=True, is_part_3=True)
                return R_all
            R_all = ceil(R_all * 1.2)  # Increment R_all

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
        return R_all

    def calc_P_X_T(self, T: int) -> tuple[Any, float, Any]:
        # Use π(T) as pi_R in Part 2 to get P_succ_single
        reconstructor = ReconstructingCombinatorialSequence(n=self.n, t=self.t, eps=self.eps, R=T, m=self.m, b=self.b,
                                                            method=self.method)
        P_succ_single, pi_T = reconstructor.calculate_probability()

        # Calculate P(X_T >= a)
        P_X_T = binom.sf(self.a - 1, self.l, P_succ_single)

        return P_X_T, P_succ_single, pi_T

    def calc_P_E(self, T, R_all):
        try:
            P_E = min(1, ((R_all + 1) ** self.l) * 2 ** (-R_all * self.epsilon(T, R_all)))
        except OverflowError:
            # Handle the overflow, for example, set result to a default value
            P_E = 1

        return P_E

    def simu_P_E(self, T, R_all):
        N = ceil(100 * 1 / self.delta)

        # Define a multinomial variable
        X = scp.stats.multinomial(ceil(R_all), [1 / self.l for _ in range(self.l)])
        # generate N draws of X
        x = X.rvs(N)
        # Find N minimal values.
        R_min = x.min(axis=1)
        # Plot distribution of minimal values
        _ = plt.hist(R_min)  # ,bins = range(50,100))
        plt.axvline(T, 0, N / 10, color='r')
        # Calculate P(E)
        P_E = np.sum(R_min < T) / N
        plt.title(f'$P(E)=P(min(R_j)<{T}) =$ {P_E}')
        plt.xlabel('$min(R_j)$')
        plt.ylabel('freq')
        plt.close()

        return P_E


if __name__ == "__main__":
    ##########
    # Part 1 #
    ##########
    n = 3  # Total number of unique building blocks in each position
    t = 2  # Required threshold on the number of observed occurrences
    eps = 0.01
    R = 40  # Acceptable error threshold

    ##########
    # Part 2 #
    ##########
    m = 10  # Total number of letters in the sequence
    b = 8  # Number of letters required to be successfully reconstructed
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
    a = ceil(l * 0.8)  # Number of barcodes required to be successfully decoded
    # eps = 0.01
    delta = 0.1  # Acceptable error threshold
    P_E_method = 'calculation'

    decoder = MessageDecoder(n=n, t=t, eps=eps, R=R, m=m, b=b, method=method, l=l, a=a, delta=delta,
                             P_E_method=P_E_method)
    R_all = decoder.run_decoder()
