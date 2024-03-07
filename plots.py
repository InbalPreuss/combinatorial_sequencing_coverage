import math
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
    t = 4  # Required threshold on the number of observed occurrences
    eps = 0.01
    R = 110  # Acceptable error threshold
    Q = 100

    number_of_simulation_runs = 50

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
        p_simu = [part1_single_position_simulator.simulate() for _ in range(number_of_simulation_runs)]
        P_simu_values.append(p_simu)
        P_simu_medians.append(np.median(p_simu))  # Calculate and store the median

        # Plotting
    plt.figure(figsize=(10, 6))

    # Line plot for calculated probability
    plt.plot(R_values, P_calc_values, label='Calculated Probability', color='blue', marker='o', linewidth=1.5)

    # Box plot for simulated values
    plt.boxplot(P_simu_values, positions=R_values, widths=1.5)

    # Line plot connecting the medians of the box plots
    plt.plot(R_values, P_simu_medians, label='Median Simulated Probability', color='red', linestyle='-', marker='x',
             linewidth=1.5)

    plt.xlabel('Number of Reads ($R$)', fontsize=20)
    # plt.tight_layout()
    plt.subplots_adjust(bottom=0.2)  # Adjust the value as needed

    plt.ylabel('Probability', fontsize=20)
    # plt.title('Probability with Different Number of Reads ($R_{\mathrm{single}}$)', fontsize=22)
    plt.xticks(R_values, fontsize=15)  # Set x-ticks to R_values
    plt.yticks(fontsize=15)  # Set x-ticks to R_values
    plt.legend(fontsize=13)
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
    R = 140  # Acceptable error threshold
    Q = 100

    # PART 2 PARAMS
    m = 100
    b = 100
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

        part2_reconstruct_combinatorial_sequence_simulator = part2_simu.CompleteSequenceSimulator(n=n, t=t, m=m, k=k,
                                                                                                  b=b, R=R_i, Q=Q)
        p_simu = [part2_reconstruct_combinatorial_sequence_simulator.simulate() for _ in range(50)]
        P_simu_values.append(p_simu)
        P_simu_medians.append(np.median(p_simu))  # Calculate and store the median

        # Plotting
    plt.figure(figsize=(10, 6))

    # Line plot for calculated probability
    plt.plot(R_values, P_calc_values, label='Calculated Probability', color='blue', marker='o', linewidth=1.5)

    # Box plot for simulated values
    plt.boxplot(P_simu_values, positions=R_values, widths=1.5)

    # Line plot connecting the medians of the box plots
    plt.plot(R_values, P_simu_medians, label='Median Simulated Probability', color='red', linestyle='-', marker='x',
             linewidth=1.5)

    plt.xlabel('Number of Reads ($R$)', fontsize=20)
    plt.subplots_adjust(bottom=0.2)  # Adjust the value as needed

    plt.ylabel('Probability', fontsize=20)
    # plt.title('Probability with Different Number of Reads (R)', fontsize=22)
    plt.xticks(R_values, fontsize=15)  # Set x-ticks to R_values
    plt.yticks(fontsize=15)  # Set x-ticks to R_values
    plt.legend(fontsize=13)
    # Save the figure as an SVG file
    plt.savefig(
        f"{PLOTS_PATH}/prob_with_different_R_combinatorial_sequence_part2_n={n},t={t},R={R},m={m},b={b}, method={method}, xaxis={R_gap}.svg",
        format='svg')
    plt.show()


def plots_part_2():
    plot_prob_with_different_R_recover_complete_combinatorial_sequence_calc_vs_simu_part2()


def plot_prob_with_different_R_recover_message_PXT_PE_part3():
    # PART 1 PARAMS
    n = 7  # Total number of unique building blocks in each position
    t = 4  # Required threshold on the number of observed occurrences
    eps = 0.01
    R = 110  # Acceptable error threshold

    # PART 2 PARAMS
    m = 10
    b = 8
    method = "binomial"

    # PART 3
    l = 10
    a = 8
    delta = 0.1
    P_E_method = 'calculation'

    R_all = 3000

    R_gap = 1
    R_values = []
    P_X_T_values = []
    P_E_values = []
    for R_i in range(45, R, R_gap):
        R_values.append(R_i)

        part3_message_decoder = part3_calc.MessageDecoder(n=n, t=t, eps=eps, m=m, b=b, method=method, l=l, a=a,
                                                          delta=delta, P_E_method=P_E_method)
        P_X_T, _, _ = part3_message_decoder.calc_P_X_T(T=R_i)
        P_X_T_values.append(P_X_T)

        # R_all = part3_message_decoder.calc_R_all(T=R_i)
        if P_E_method == 'simulation':
            P_E = part3_message_decoder.simu_P_E(R_all=R_all, T=R_i)
        else:
            P_E = part3_message_decoder.calc_P_E(R_all=R_all, T=R_i)
        print(f'1-P_E={1 - P_E}')
        P_E_values.append(1 - P_E)

        # Plotting
    plt.figure(figsize=(10, 6))

    # Line plot for calculated probability
    plt.plot(R_values, P_X_T_values, label=r'$P_{\mathrm{X}_{\mathrm{ρ}}}$', color='blue', marker='o')
    plt.plot(R_values, P_E_values, label='1-$P(E)$', color='red', marker='o')

    # Delta
    threshold = math.sqrt(1-delta)
    plt.axhline(y=threshold, color='deeppink', linestyle='--')
    # plt.text(R_values[-1], delta_example, 'sqrt(1-δ)', color='red', verticalalignment='bottom')
    plt.text(max(R_values) + 5, threshold, r'$\sqrt{1-\delta}$', color='deeppink', verticalalignment='bottom',
             fontsize=12)

    # Find indices where y-values are >= delta_example and add vertical lines
    indices_P_E = [i for i, y in enumerate(P_E_values) if y >= threshold]
    indices_P_X_T = [i for i, y in enumerate(P_X_T_values) if y >= threshold]
    if indices_P_X_T and indices_P_E:
        first_idx, last_idx = indices_P_X_T[0], indices_P_E[-1]
        plt.axvline(x=R_values[first_idx], color='deeppink', linestyle='--')
        plt.axvline(x=R_values[last_idx], color='deeppink', linestyle='--')

    # Add horizontal line at y = 0.8
    delta_example = math.sqrt(1-0.2)
    plt.axhline(y=delta_example, color='green', linestyle='--')
    # plt.text(R_values[-1], delta_example, 'sqrt(1-δ)', color='red', verticalalignment='bottom')
    plt.text(max(R_values) + 5, delta_example, r'$\sqrt{1-\delta}$', color='green', verticalalignment='bottom',
             fontsize=12)

    # Find indices where y-values are >= delta_example and add vertical lines
    indices_P_E = [i for i, y in enumerate(P_E_values) if y >= delta_example]
    indices_P_X_T = [i for i, y in enumerate(P_X_T_values) if y >= delta_example]
    if indices_P_X_T and indices_P_E:
        first_idx, last_idx = indices_P_X_T[0], indices_P_E[-1]
        plt.axvline(x=R_values[first_idx], color='green', linestyle='--')
        plt.axvline(x=R_values[last_idx], color='green', linestyle='--')

    plt.xlabel('Number of Reads (ρ)', fontsize=20)
    plt.ylabel('Decoding Probability', fontsize=20)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    # plt.title('Probability with Different Number of Reads (T)')
    plt.xticks(range(min(R_values), max(R_values)+1, 5))  # Set x-ticks to R_values
    plt.legend(fontsize=13)
    # plt.grid(True)
    # Save the figure as an SVG file
    plt.savefig(
        f"{PLOTS_PATH}/prob_with_different_R_recover_message_PXT_PE_part3_n={n},t={t},R={R},m={m},b={b},l={l},a={a},delta={delta},P_E_method={P_E_method}, xaxis={R_gap}.svg",
        format='svg')
    plt.show()


def plots_part_3():
    plot_prob_with_different_R_recover_message_PXT_PE_part3()


if __name__ == '__main__':
    if not os.path.exists(PLOTS_PATH):
        os.makedirs(PLOTS_PATH)

    plots_part_1()
    plots_part_2()
    plots_part_3()
