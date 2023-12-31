import numpy as np
from multiprocessing import Pool

from simulation.part2_reconstructing_a_complete_combinatorial_sequence_simulation import CompleteSequenceSimulator

class CompleteMessageSimulator:
    def __init__(self, n, t, m, k, b, l, a, R_all, Q):
        self.n = n
        self.t = t
        self.m = m
        self.k = k
        self.b = b
        self.l = l
        self.a = a
        self.R_all = R_all
        self.Q = Q

    def simulate(self):
        success_count = 0
        for _ in range(self.Q):
            reads_distribution = np.random.multinomial(self.R_all, [1/self.l]*self.l)
            decoded_sequences = 0
            for R_j in reads_distribution:
                if CompleteSequenceSimulator(self.n, self.t, self.m, self.k, self.b, R_j, 1).simulate():
                    decoded_sequences += 1
            if decoded_sequences >= self.a:
                success_count += 1
        return success_count / self.Q

    def simulate_single_sequence(self, n, t, m, k, b, R_j):
        sequence_simulator = CompleteSequenceSimulator(n, t, m, k, b, R_j, 1)
        return sequence_simulator.simulate()

    def simulate_parallel(self):
        success_count = 0
        with Pool() as pool:
            for _ in range(self.Q):
                reads_distribution = np.random.multinomial(self.R_all, [1 / self.l] * self.l)
                results = pool.starmap(self.simulate_single_sequence,
                                       [(self.n, self.t, self.m, self.k, self.b, R_j) for R_j in reads_distribution])
                decoded_sequences = sum(results)
                if decoded_sequences >= self.a:
                    success_count += 1
        return success_count / self.Q





