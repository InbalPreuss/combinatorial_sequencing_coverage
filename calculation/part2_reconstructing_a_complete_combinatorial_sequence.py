from typing import Tuple, Any

import numpy as np
from scipy.stats import binom, norm

from calculation.part1_reconstructing_a_single_combinatorial_position import ReconstructingSingleCombinatorialPosition


class ReconstructingCombinatorialSequence:

    def __init__(self,n: int, t: int, eps: float, R: int, m: int, b: int, method: str):
        """
        Initialize the parameters for reconstructing the combinatorial sequence.

        :param m: Total number of letters in the combinatorial sequence.
        :param b: Number of letters required to be successfully reconstructed.
        :param pi_R: Probability of successfully reconstructing a single letter.
        """
        self.n = n
        self.t = t
        self.eps = eps
        self.R = R
        self.m = m
        self.b = b
        self.pi_R = 0.0
        self.method = method

    def calculate_probability(self) -> tuple[float, Any]:
        self.pi_R = self.run_part1_for_prob()
        """
        Calculate the probability of successfully decoding the sequence.

        :return: Probability of successful decoding.
        """
        if self.method == 'binomial':
            return self._binomial_method(), self.pi_R
        elif self.method == 'normal':
            return self._normal_approximation(), self.pi_R
        else:
            raise ValueError("Invalid method. Choose 'binomial' or 'normal'.")

    def _binomial_method(self) -> float:
        """ Use the binomial distribution formula for calculation. """
        # prob = sum(binom.pmf(k, self.m, self.pi_R) for k in range(self.b, self.m + 1))
        # TODO: need to call part1 with the params so we can get the pi_R form the calculation
        prob = binom.sf(self.b - 1, self.m, self.pi_R)
        return prob

    def _normal_approximation(self) -> float:
        """ Use the normal approximation for calculation. """
        mean = self.m * self.pi_R
        std_dev = np.sqrt(self.m * self.pi_R * (1 - self.pi_R))
        prob = 1 - norm.cdf(self.b - 1, mean, std_dev)
        return prob

    def run_part1_for_prob(self):
        part1_reconstructing_single_combinatorial_position = ReconstructingSingleCombinatorialPosition(self.n, self.t, self.eps, self.R)
        _, p, _ = part1_reconstructing_single_combinatorial_position.calc_results()

        return p

if __name__ == '__main__':
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
    m = 10  # Total number of letters in the sequence
    b = 9  # Number of letters required to be successfully reconstructed
    # pi_R = 0.8  # Probability of successfully reconstructing a single letter
    method = "normal"
    reconstructor = ReconstructingCombinatorialSequence(n=n, t=t, eps=eps, R=R, m=m, b=b, method=method)
    P_succ_single, pi_R = reconstructor.calculate_probability()
