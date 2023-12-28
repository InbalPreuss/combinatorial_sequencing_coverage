import numpy as np
import pandas as pd
# from scipy.special import binom
import matplotlib.pyplot as plt
from itertools import combinations
import sympy
from sympy.abc import pi

# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    import numpy as np
    import pandas as pd
    # from scipy.special import binom
    import matplotlib.pyplot as plt
    from itertools import combinations
    import sympy
    from sympy.abc import pi


    class General_Coupon_Collector():

        def set_params(self, n, k, t, eps):
            """

            :param n: Number of unique coupons.
            :param k: Number of unique coupons to collect.
            :param t: Number of copies of unique coupon to collect.
            """
            self.n = n
            self.k = k
            self.t = t
            self.eps = eps
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
                    p = (1-self.eps) * (state[j] / self.n)
                    new_state = state.copy()
                    new_state[j] -= 1
                    new_state[j + 1] += 1
                    new_state_idx = self.states.index(new_state)
                    self.transition_matrix[state_idx, new_state_idx] = p
                p = (1-self.eps) * (state[self.t] / self.n) + self.eps
                self.transition_matrix[state_idx, state_idx] = p

        def _calc_init_vec(self):
            self.init_vec = np.zeros(self.n_state)
            self.init_vec[0] = 1

        def calc_results(self, s):
            """
            :param s: Number of rounds for collecting the copupon.
            :return:
            """
            pi = "\u03c0"
            res_dist = self.init_vec @ np.linalg.matrix_power(self.transition_matrix, s)
            res_prob = res_dist[np.array(self.states)[:, self.t] >= self.k].sum()
            arr = np.column_stack((np.array([np.array(x) for x in self.states]), np.array(res_dist)))
            dict_for_plot = {tuple(x[:t + 1].astype(int)): x[-1] for x in \
                             pd.DataFrame(arr).sort_values([x for x in range(t, 0, -1)]).values}
            fig = plt.figure(figsize=(15, 8))
            plt.bar([str(x) for x in dict_for_plot.keys()], dict_for_plot.values())
            labels = [x if i % 1 == 0 else '' for i, x in enumerate(dict_for_plot.keys())]
            plt.ylim(0, 1)
            plt.ylabel('Probability', fontsize=20)
            plt.xlabel('States', fontsize=20)
            plt.xticks(fontsize=15)
            plt.yticks(fontsize=15)
            plt.title(
                'Coupon Collector - ${}(R={})=$P(T(n={},t={})) $\leq {})=${:.3f}'.format(pi, s, self.n,
                                                                                              self.t, s, res_prob))

            return dict_for_plot, res_prob, {tuple(k): v for k, v in zip(self.states, res_dist)}, [
                [str(x) for x in dict_for_plot.keys()], dict_for_plot.values()]


    from matplotlib.animation import FuncAnimation
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 1)
    fig.set_size_inches(5, 5)

    n = 5
    k = n
    t = 1
    eps = 0
    gcc = General_Coupon_Collector()
    gcc.set_params(n, k, t, eps)
    s = 30

    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.animation import FuncAnimation

    n_frames = s


    def update(frame):
        dict_for_plot, _, _, plot_data = gcc.calc_results(frame)
        ax.clear()
        ax.bar(plot_data[0], plot_data[1])
        ax.set_ylim(ymin=0, ymax=1)
        ax.set_title('Coupon Collector - s={}'.format(frame))


    gcc = General_Coupon_Collector()
    gcc.set_params(n=n, k=n, t=t)

    ani = FuncAnimation(fig, update, frames=range(1, n_frames + 1))
    ani.save('animation.gif', writer='ffmpeg')

    for frame in range(1, s + 1):
        fig, ax = plt.subplots(1, 1)
        fig.set_size_inches(5, 5)

        dict_for_plot, _, _, plot_data = gcc.calc_results(frame)
        ax.bar(plot_data[0], plot_data[1])
        ax.set_ylim(ymin=0, ymax=1)
        ax.set_title('Coupon Collector - s={}'.format(frame))

        # Save the figure
        plt.savefig(f'plot_{frame}.svg', format='svg')
        plt.close(fig)  # Close the figure to free memory

