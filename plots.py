import os

import numpy as np
from matplotlib import pyplot as plt

# PART 1
import calculation.part1_reconstructing_a_single_combinatorial_position as part1_calc
import simulation.part1_reconstructing_a_single_combinatorial_position_simulation as part1_simu

# PART 2
import calculation.part2_reconstructing_a_complete_combinatorial_sequence as part2_calc
import simulation.part2_reconstructing_a_complete_combinatorial_sequence_simulation as part2_simu

# PART 3
import calculation.part3_reconstructing_a_complete_message as part3_calc
import simulation.part3_reconstructing_a_complete_message_simulation as part3_simu

PLOTS_PATH = "plots"


def plot_prob_with_different_R_recover_single_letter_calc_vs_simu_part1():
    n = 7  # Total number of unique building blocks in each position
    t = 1  # Required threshold on the number of observed occurrences
    eps = 0.01
    R = 100  # Acceptable error threshold
    Q = 100

    R_values = []
    P_calc_values = []
    P_simu_values = []
    P_simu_medians = []
    for R_i in range(10, R, 10):
        R_values.append(R_i)

        part1_reconstructing_single_combinatorial_position = part1_calc.ReconstructingSingleCombinatorialPosition(n=n,
                                                                                                                  t=t,
                                                                                                                  eps=eps,
                                                                                                                  R=R_i)
        _, p_calc, _ = part1_reconstructing_single_combinatorial_position.calc_results()
        P_calc_values.append(p_calc)

        part1_single_position_simulator = part1_simu.SinglePositionSimulator(n=n, t=t, R=R_i, Q=Q)
        p_simu = [part1_single_position_simulator.simulate() for _ in range(50)]
        P_simu_values.append(p_simu)
        P_simu_medians.append(np.median(p_simu))  # Calculate and store the median

        # Plotting
    plt.figure(figsize=(10, 6))

    # Line plot for calculated probability
    plt.plot(R_values, P_calc_values, label='Calculated Probability', color='blue', marker='o')

    # Box plot for simulated values
    plt.boxplot(P_simu_values, positions=R_values, widths=1.0)

    # Line plot connecting the medians of the box plots
    plt.plot(R_values, P_simu_medians, label='Median Simulated Probability', color='red', linestyle='-', marker='x')

    plt.xlabel('Number of Reads (R)')
    plt.ylabel('Probability')
    plt.title('Probability with Different Number of Reads (R)')
    plt.xticks(R_values)  # Set x-ticks to R_values
    plt.legend()
    plt.grid(True)
    # Save the figure as an SVG file
    plt.savefig(f"{PLOTS_PATH}/prob_with_different_R_recover_single_letter_part1_n={n},t={t},R={R}.svg", format='svg')
    plt.show()


def plots_part_1():
    plot_prob_with_different_R_recover_single_letter_calc_vs_simu_part1()


def plot_prob_with_different_R_recover_complete_combinatorial_sequence_calc_vs_simu_part2():
    # PART 1 PARAMS
    n = 7  # Total number of unique building blocks in each position
    t = 4  # Required threshold on the number of observed occurrences
    eps = 0.01
    R = 130  # Acceptable error threshold
    Q = 100

    # PART 2 PARAMS
    m = 100
    b = 95
    method = "binomial"
    k = n

    R_gap = 10
    R_values = []
    P_calc_values = []
    P_simu_values = []
    P_simu_medians = []
    for R_i in range(10, R, R_gap):
        R_values.append(R_i)

        part2_reconstruct_combinatorial_sequence = part2_calc.ReconstructingCombinatorialSequence(n=n, t=t, eps=eps,
                                                                                                  R=R_i,
                                                                                                  m=m, b=b,
                                                                                                  method=method)
        p_calc, _ = part2_reconstruct_combinatorial_sequence.calculate_probability()
        P_calc_values.append(p_calc)

        part2_reconstruct_combinatorial_sequence_simulator = part2_simu.CompleteSequenceSimulator(n=n, t=t, m=m, k=k, b=b, R=R_i, Q=Q)
        p_simu = [part2_reconstruct_combinatorial_sequence_simulator.simulate() for _ in range(50)]
        P_simu_values.append(p_simu)
        P_simu_medians.append(np.median(p_simu))  # Calculate and store the median

        # Plotting
    plt.figure(figsize=(10, 6))

    # Line plot for calculated probability
    plt.plot(R_values, P_calc_values, label='Calculated Probability', color='blue', marker='o')

    # Box plot for simulated values
    plt.boxplot(P_simu_values, positions=R_values, widths=1.0)

    # Line plot connecting the medians of the box plots
    plt.plot(R_values, P_simu_medians, label='Median Simulated Probability', color='red', linestyle='-', marker='x')

    plt.xlabel('Number of Reads (R)')
    plt.ylabel('Probability')
    plt.title('Probability with Different Number of Reads (R)')
    plt.xticks(R_values)  # Set x-ticks to R_values
    plt.legend()
    plt.grid(True)
    # Save the figure as an SVG file
    plt.savefig(f"{PLOTS_PATH}/prob_with_different_R_combinatorial_sequence_part2_n={n},t={t},R={R},m={m},b={b},xaxis={R_gap}.svg", format='svg')
    plt.show()


def plots_part_2():
    plot_prob_with_different_R_recover_complete_combinatorial_sequence_calc_vs_simu_part2()


def plot_prob_with_different_R_recover_message_PXT_PE_part3():
    # PART 1 PARAMS
    n = 7  # Total number of unique building blocks in each position
    t = 4  # Required threshold on the number of observed occurrences
    eps = 0.01
    R = 96  # Acceptable error threshold

    # PART 2 PARAMS
    m = 100
    b = 80
    method = "binomial"

    # PART 3
    l = 50
    a = 40
    delta = 0.1

    R_gap = 1
    R_values = []
    P_X_T_values = []
    P_E_values = []
    for R_i in range(55, R, R_gap):
        R_values.append(R_i)

        part3_message_decoder = part3_calc.MessageDecoder(n=n, t=t, eps=eps, R=R_i,m=m, b=b, method=method, l=l, a=a, delta=delta)
        P_X_T, _, _ = part3_message_decoder.calc_P_X_T(T=R_i)
        P_X_T_values.append(P_X_T)

        R_all = 55000
        # R_all = part3_message_decoder.calc_R_all(T=R_i)
        P_E = part3_message_decoder.calc_P_E(R_all=R_all, T=R_i)
        print(f'1-P_E={1-P_E}')
        P_E_values.append(1-P_E)

        # Plotting
    plt.figure(figsize=(10, 6))

    # Line plot for calculated probability
    plt.plot(R_values, P_X_T_values, label='P_X_T', color='blue', marker='o')
    plt.plot(R_values, P_E_values, label='PE', color='red', marker='o')



    plt.xlabel('Number of Reads (R)')
    plt.ylabel('Probability')
    plt.title('Probability with Different Number of Reads (R)')
    plt.xticks(R_values)  # Set x-ticks to R_values
    plt.legend()
    plt.grid(True)
    # Save the figure as an SVG file
    plt.savefig(
        f"{PLOTS_PATH}/prob_with_different_R_recover_message_PXT_PE_part3_n={n},t={t},R={R},m={m},b={b},l={l},a={a},delta={delta},xaxis={R_gap}.svg",
        format='svg')
    plt.show()


def plots_part_3():
    plot_prob_with_different_R_recover_message_PXT_PE_part3()


if __name__ == '__main__':
    if not os.path.exists(PLOTS_PATH):
        os.makedirs(PLOTS_PATH)

    # plots_part_1()
    # plots_part_2()
    plots_part_3()
