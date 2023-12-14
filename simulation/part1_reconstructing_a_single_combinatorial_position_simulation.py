import numpy as np
from matplotlib import pyplot as plt


class SinglePositionSimulator:
    def __init__(self, n, t, R, Q):
        self.n = n
        self.t = t
        self.R = R
        self.Q = Q

    def simulate(self):
        success_count = 0
        for _ in range(self.Q):
            reads = np.random.choice(range(1, self.n + 1), self.R, replace=True)
            if all(np.count_nonzero(reads == i) >= self.t for i in range(1, self.n + 1)):
                success_count += 1
        return success_count / self.Q

    def visualize_values_and_prob(self, values, x_label, P_values, y_label, title):
        plt.plot(values, P_values, marker='o')
        plt.xlabel(x_label)
        plt.ylabel(y_label)
        plt.title(f'Finding {title}')
        plt.show()
