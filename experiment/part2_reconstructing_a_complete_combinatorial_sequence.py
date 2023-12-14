import numpy as np
from scipy.stats import binom, norm


class ReconstructingCombinatorialSequence:

    def __init__(self, m: int, b: int, pi_R: float, method: str):
        """
        Initialize the parameters for reconstructing the combinatorial sequence.

        :param m: Total number of letters in the combinatorial sequence.
        :param b: Number of letters required to be successfully reconstructed.
        :param pi_R: Probability of successfully reconstructing a single letter.
        """
        self.m = m
        self.b = b
        self.pi_R = pi_R
        self.method = method

    def calculate_probability(self) -> float:
        """
        Calculate the probability of successfully decoding the sequence.

        :return: Probability of successful decoding.
        """
        if self.method == 'binomial':
            return self._binomial_method()
        elif self.method == 'normal':
            return self._normal_approximation()
        else:
            raise ValueError("Invalid method. Choose 'binomial' or 'normal'.")

    def _binomial_method(self) -> float:
        """ Use the binomial distribution formula for calculation. """
        prob = sum(binom.pmf(k, self.m, self.pi_R) for k in range(self.b, self.m + 1))
        return prob

    def _normal_approximation(self) -> float:
        """ Use the normal approximation for calculation. """
        mean = self.m * self.pi_R
        std_dev = np.sqrt(self.m * self.pi_R * (1 - self.pi_R))
        prob = 1 - norm.cdf(self.b - 1, mean, std_dev)
        return prob


# Example usage
def run_part2_reconstructing_combinatorial_sequence(m: int, b: int, pi_R: float, method: str) -> float:
    reconstructor = ReconstructingCombinatorialSequence(m, b, pi_R, method=method)
    P_succ_single = reconstructor.calculate_probability()

    print(f"Probability of successful decoding: {P_succ_single}")

    return P_succ_single


if __name__ == '__main__':
    m = 10  # Total number of letters in the sequence
    b = 7  # Number of letters required to be successfully reconstructed
    pi_R = 0.8  # Probability of successfully reconstructing a single letter
    method = "normal"

    P_succ_single = run_part2_reconstructing_combinatorial_sequence(m=m,b=b,pi_R=pi_R,method=method)
