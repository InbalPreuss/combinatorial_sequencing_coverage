import os
from itertools import product
import numpy as np
import pandas as pd
from tqdm import tqdm

from part1_reconstructing_a_single_combinatorial_position_simulation import SinglePositionSimulator
from part2_reconstructing_a_complete_combinatorial_sequence_simulation import CompleteSequenceSimulator
from part3_reconstructing_a_complete_message_simulation import CompleteMessageSimulator

# Define parameters
n_values = [3, 5, 7]
t_values = [1, 2, 3, 4]
R_values = list(range(10, 300, 10))
Q_values = [300]

k_values = [3, 5, 7]
m_values = [5, 10, 20, 40, 60, 80, 100, 150, 200]
b_fractions = [1, 0.975, 0.95, 0.9, 0.85, 0.8, 0.5]

l_values = [10, 20, 30, 40, 50, 100, 150, 200]
a_fractions = [1, 0.975, 0.95, 0.9, 0.85, 0.8, 0.5]

R_all_range = list(range(5, 100, 5))
R_all_values = [np.floor(np.array(R_all_range) * l).astype(int) for l in l_values]

if __name__ == '__main__':
    # Directory path
    output_dir = 'output/'

    # Create the directory
    os.makedirs(output_dir, exist_ok=True)

    # ###################################
    # ###### Run Part 1 Simulation ######
    # ###################################
    # #TODO: delete the examples
    # # n_values, t_values, R_values, Q_values = [1], [1], [1, 2, 3], [1]
    #
    # # Calculate total iterations for the progress bar
    # total_iterations = len(n_values) * len(t_values) * len(R_values) * len(Q_values)
    #
    # # CSV preparation
    # part1_csv_path = os.path.join(output_dir, 'part1_SinglePositionSimulatorResults.csv')
    # part1_headers = ['n', 't', 'R', 'Q', 'Result']
    #
    # # Open the CSV file and write headers
    # with open(part1_csv_path, 'w', newline='') as file:
    #     writer = pd.DataFrame(columns=part1_headers)
    #     writer.to_csv(file, index=False)
    #     for n, t, R, Q in tqdm(product(n_values, t_values, R_values, Q_values), total=total_iterations):
    #         part1_simulator = SinglePositionSimulator(n, t, R, Q)
    #         part1_result = part1_simulator.simulate()
    #
    #         df_row = pd.DataFrame([[n, t, R, Q, part1_result]], columns=part1_headers)
    #         df_row.to_csv(file, mode='a', header=False, index=False)
    #

    # ###################################
    # ###### Run Part 2 Simulation ######
    # ####################################
    # # Calculate total iterations for the progress bar
    # total_iterations = len(n_values) * len(t_values) * len(R_values) * len(Q_values) * len(m_values) * len(
    #     b_fractions) * len(k_values)
    #
    # # log file
    # progress_log_path = os.path.join(output_dir, 'progress2.log')
    # # Write the total number of iterations to the log file
    # with open(progress_log_path, 'a') as progress_file:
    #     progress_file.write(f"Total Iterations: {total_iterations}\n")
    #
    # # CSV preparation
    # part2_csv_path = os.path.join(output_dir, 'part2_CompleteSequenceSimulator.csv')
    # part2_headers = ['n', 't', 'R', 'm', 'b', 'Q', 'k', 'Result']
    #
    # with open(part2_csv_path, 'w', newline='') as file:
    #     writer = pd.DataFrame(columns=part2_headers)
    #     writer.to_csv(file, index=False)
    #     count_iter = 0
    #     for n, t, R, m, b_fraction, Q, k in tqdm(product(n_values, t_values, R_values, m_values, b_fractions, Q_values, k_values), total=total_iterations):
    #         count_iter += 1
    #         if count_iter % 100 == 0:
    #             with open(progress_log_path, 'a') as progress_file:
    #                 progress_file.write(f"Iteration number: {count_iter}\n")
    #
    #         b = int(np.floor(b_fraction * m))
    #
    #         if n < k:
    #             continue
    #         part2_simulator = CompleteSequenceSimulator(n, t, m, k, b, R, Q)
    #         part2_result = part2_simulator.simulate_parallel()
    #         # part2_result = part2_simulator.simulate()
    #
    #         df_row = pd.DataFrame([[n, t, R, m, b, Q, k, part2_result]], columns=part2_headers)
    #         df_row.to_csv(file, mode='a', header=False, index=False)
    #         file.flush()
    #
    #         # print(f"Part 2 Result: {part2_result}")
    #         # print(f"n={n}, t={t}, R={R}, Q={Q}, b={b}, k={k}")
    #     with open(progress_log_path, 'w') as progress_file:
    #         progress_file.write(f"End of Iterations part 2\n")

    ###################################
    ###### Run Part 3 Simulation ######
    ###################################

    # Calculate total iterations for the progress bar
    total_iterations = len(n_values) * len(t_values) * len(R_values) * len(Q_values) * len(m_values) * len(
        b_fractions) * len(k_values) * len(l_values) * len(a_fractions) * len(R_all_range)

    print(total_iterations)
    # CSV preparation
    part3_csv_path = os.path.join(output_dir, 'part3_CompleteMessageSimulator.csv')
    part3_headers = ['n', 't', 'm', 'b', 'Q', 'k', 'l', 'a', 'R_all', 'Result']
    # Open the CSV file and write headers
    with open(part3_csv_path, 'w', newline='') as file:
        writer = pd.DataFrame(columns=part3_headers)
        writer.to_csv(file, index=False)
        for n, t, m, b_fraction, Q, k, l, a_fraction, R_all_value in tqdm(
                product(n_values, t_values, m_values, b_fractions, Q_values, k_values, l_values, a_fractions,
                        R_all_range), total=total_iterations):
            a = int(np.floor(a_fraction * l))
            R_all = int(np.floor(R_all_value * l))
            b = int(np.floor(b_fraction * m))

            if n < k:
                continue
            part3_simulator = CompleteMessageSimulator(n, t, m, k, b, l, a, R_all, Q)
            part3_result = part3_simulator.simulate_parallel()

            df_row = pd.DataFrame([[n, t, m, b, Q, k, l, a, R_all, part3_result]], columns=part3_headers)
            df_row.to_csv(file, mode='a', header=False, index=False)
            file.flush()

            # print(f"Part 3 Result: {part3_result}")
            # print(f"n={n}, t={t}, R_all={R_all}, Q={Q}, b={b}, k={k}, m={m}, a={a},l={l}")
