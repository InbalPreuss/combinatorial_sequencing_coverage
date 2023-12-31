import scipy as scp
import numpy as np
import matplotlib.pyplot as plt
import os

def make_dir(path):
    if not os.path.exists(path):
        os.makedirs(path)

if __name__ == '__main__':
    R_all = 4500
    l = 50
    N = 100000
    T = 65
    X = scp.stats.multinomial(R_all, [1/l for _ in range(l)])
    x = X.rvs(N)
    R_min = x.min(axis=1)

    plt.figure(figsize=(10,6))
    plt.subplots_adjust(bottom=0.2)
    counts, bins, patches = plt.hist(R_min, bins=range(min(R_min), max(R_min) + 1, 1), edgecolor='black')

    # Change color of bins
    for patch, leftside in zip(patches, bins):
        if leftside < T:
            patch.set_facecolor('gray')
        else:
            patch.set_facecolor('blue')

    # Marking P(E) and 1-P(E)
    plt.axvline(T, color='r', label='P(E) Threshold at ρ')

    P_E = np.sum(R_min < T) / N
    # Adding a text box for P(E)
    textstr = f'P(E) = {P_E:.4f}'
    props = dict(boxstyle='round', facecolor='gray', alpha=0.5)
    plt.text(0.05, 0.95, textstr, transform=plt.gca().transAxes, fontsize=14,
             verticalalignment='top', bbox=props)

    # plt.title(f'Error Probability P(E) = {P_E}')
    plt.xlabel('Minimum Reads ($R_{\mathrm{min}}$)', fontsize=20)
    plt.ylabel('Frequency', fontsize=20)
    plt.legend()
    plt.yticks(fontsize=15)
    plt.xticks(fontsize=15)
    plt.xlim(45,95)

    # Update with your PLOTS_PATH
    PLOTS_PATH = 'plots/'
    make_dir(PLOTS_PATH)
    plt.savefig(f"{PLOTS_PATH}/Multinomial_bound_R_all={R_all}_T={T}_l={l}_N={N}.svg", format='svg')
    plt.show()
