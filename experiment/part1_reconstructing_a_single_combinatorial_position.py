import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from itertools import combinations


class ReconstructingSingleCombinatorialPosition:

    def __init__(self, n, t, eps, R):
        """
        Set the parameters for the reconstruction simulation.

        :param n: Number of unique coupons.
        :param t: Number of copies of unique coupon to collect.
        :param eps: Error probability.
        """
        self.n = n
        self.t = t
        self.eps = eps
        self.R = R
        self._calc_state()
        self._calc_transition_matrix()
        self._calc_init_vec()

    def _calc_state(self):
        self.states = []
        for c in list(combinations(range(self.n + self.t), self.t))[::-1]:
            self.states.append([b - a - 1 for a, b in zip((-1,) + c, c + (self.n + self.t,))])
        self.n_state = len(self.states)

    def _calc_transition_matrix(self):
        self.transition_matrix = np.zeros((self.n_state, self.n_state))
        for state_idx, state in enumerate(self.states):
            for j in range(self.t):
                if state[j] == 0:
                    continue
                p = (1 - self.eps) * (state[j] / self.n)
                new_state = state.copy()
                new_state[j] -= 1
                new_state[j + 1] += 1
                new_state_idx = self.states.index(new_state)
                self.transition_matrix[state_idx, new_state_idx] = p
            p = (1 - self.eps) * (state[self.t] / self.n) + self.eps
            self.transition_matrix[state_idx, state_idx] = p

    def _calc_init_vec(self):
        self.init_vec = np.zeros(self.n_state)
        self.init_vec[0] = 1

    def calc_results(self):
        """
        Calculate the results of the reconstruction simulation.


        :param R: Number of rounds for collecting the copupon.
        :return:
        """
        pi = "\u03c0"
        res_dist = self.init_vec @ np.linalg.matrix_power(self.transition_matrix, self.R)
        res_prob = res_dist[np.array(self.states)[:, self.t] >= self.n].sum()
        arr = np.column_stack((np.array([np.array(x) for x in self.states]), np.array(res_dist)))
        dict_for_plot = {tuple(x[:self.t + 1].astype(int)): x[-1] for x in \
                         pd.DataFrame(arr).sort_values([x for x in range(self.t, 0, -1)]).values}
        fig = plt.figure(figsize=(15, 8))
        plt.bar([str(x) for x in dict_for_plot.keys()], dict_for_plot.values())
        labels = [x if i % 1 == 0 else '' for i, x in enumerate(dict_for_plot.keys())]
        plt.xticks(range(len(self.states)), labels, rotation='vertical')
        plt.ylabel('Probability')
        plt.xlabel('State')
        plt.title(
            'Coupon Collector - ${}(s={})=$P(T(n={},t={})) $\leq {})=${:.3f}'.format(pi, self.R, self.n, self.t,
                                                                                     self.R, res_prob))
        # Red line
        for i, state in enumerate(self.states):
            if state[self.t] >= self.n:
                index = list(dict_for_plot).index(tuple(state))
                break
        plt.vlines(index, 0, max(res_dist), 'r')

        plt.show()
        return dict_for_plot, res_prob, {tuple(k): v for k, v in zip(self.states, res_dist)}


def run_part1_reconstructing_single_combinatorial_position(n, t, eps, R):
    gcc = ReconstructingSingleCombinatorialPosition(n, t, eps, R)
    rr, p, r = gcc.calc_results()
    return rr, p, r

if __name__ == '__main__':
    n = 10
    t = 3
    eps = 0.1
    R = 40

    n = 5  # Total number of unique building blocks in each position
    t = 3  # Required threshold on the number of observed occurrences
    eps = 0.01
    R = 40  # Acceptable error threshold
    run_part1_reconstructing_single_combinatorial_position(n, t, eps, R)
