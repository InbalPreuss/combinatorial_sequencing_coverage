import os

import numpy as np
from multiprocessing import Pool

from simulation.part1_reconstructing_a_single_combinatorial_position_simulation import SinglePositionSimulator

class CompleteSequenceSimulator:
    def __init__(self, n, t, m, k, b, R, Q):
        self.n = n
        self.t = t
        self.m = m
        self.k = k
        self.b = b
        self.R = R
        self.Q = Q

    # not parallel
    def simulate(self):
        success_count = 0
        for _ in range(self.Q):
            # sequence = [np.random.choice(range(1, self.n + 1), self.k, replace=False) for _ in range(self.m)]
            # decoded_count = sum(SinglePositionSimulator(self.n, self.t, self.R, self.Q).simulate() for _ in sequence)
            decoded_count = sum(SinglePositionSimulator(self.k, self.t, self.R, 1).simulate() for _ in range(self.m))
            if decoded_count >= self.b:
                success_count += 1
        return success_count / self.Q

    def single_run_simulation(self, n, t, R, m, k, b):
        single_simulator = SinglePositionSimulator(n, t, R, 1)
        # sequence = [np.random.choice(range(1, n + 1), k, replace=False) for _ in range(m)]
        # decoded_count = sum(single_simulator.simulate() for _ in sequence)
        decoded_count = sum(single_simulator.simulate() for _ in range(self.m))
        return int(decoded_count >= b)

    def simulate_parallel(self):
        # Use one fewer core than the total available
        num_cores = max(os.cpu_count() - 1, 1)  # Ensure at least one core is used
        with Pool(processes=num_cores) as pool:
            results = pool.starmap(self.single_run_simulation,
                                   [(self.n, self.t, self.R, self.m, self.k, self.b) for _ in range(self.Q)])
        success_count = sum(results)
        return success_count / self.Q

