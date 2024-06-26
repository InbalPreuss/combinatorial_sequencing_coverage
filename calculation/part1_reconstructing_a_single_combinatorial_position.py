import re

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from itertools import combinations

DIR_PLOT_PATH = 'plots/'
DIR_PLOT_PATH_PART_1 = f'{DIR_PLOT_PATH}calculation_part1_plots/'

import utils as uts
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
        uts.make_dir(DIR_PLOT_PATH)
        uts.make_dir(DIR_PLOT_PATH_PART_1)

        pi = "\u03c0"
        res_dist = self.init_vec @ np.linalg.matrix_power(self.transition_matrix, self.R)
        res_prob = min(res_dist[np.array(self.states)[:, self.t] >= self.n].sum(), 0.99)
        arr = np.column_stack((np.array([np.array(x) for x in self.states]), np.array(res_dist)))
        dict_for_plot = {tuple(x[:self.t + 1].astype(int)): x[-1] for x in \
                         pd.DataFrame(arr).sort_values([x for x in range(self.t, 0, -1)]).values}
        plt.figure(figsize=(10, 6))
        plt.bar([str(x) for x in dict_for_plot.keys()], dict_for_plot.values())
        labels = [x if i % 1 == 0 else '' for i, x in enumerate(dict_for_plot.keys())]
        plt.xticks(range(len(self.states)), labels, rotation='vertical', fontsize=15)
        plt.yticks(fontsize=15)
        plt.ylabel('Probability', fontsize=20)
        plt.xlabel('States ($S=(v(0),v(1),v(2))$)', fontsize=20)
        plt.ylim(0,1)
        plt.subplots_adjust(bottom=0.23)
        title = 'Coupon Collector - ${}(R={})=$P(T(n={},t={})) $\leq {})=${:.3f}'.format(pi, self.R, self.n, self.t, self.R, res_prob)
        # plt.title(title)
        # Red line
        for i, state in enumerate(self.states):
            if state[self.t] >= self.n:
                index = list(dict_for_plot).index(tuple(state))
                break
        plt.vlines(index, 0, max(res_dist), 'r')

        safe_title = f'Coupon Collector_pi(R={self.R})=P(T(n={self.n},t={self.t}))smallerthen{self.R}={res_prob})'
        plt.savefig(f"{DIR_PLOT_PATH_PART_1}/{safe_title}.svg",
                    format='svg')

        # plt.show()
        plt.close()
        return dict_for_plot, res_prob, {tuple(k): v for k, v in zip(self.states, res_dist)}

if __name__ == '__main__':
    n = 3  # Total number of unique building blocks in each position
    t = 2  # Required threshold on the number of observed occurrences
    eps = 0.01
    R = 5  # Acceptable error threshold
    part1_reconstructing_single_combinatorial_position = ReconstructingSingleCombinatorialPosition(n, t, eps, R)
    rr, p, r = part1_reconstructing_single_combinatorial_position.calc_results()
