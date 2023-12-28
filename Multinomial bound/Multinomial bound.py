# import scipy as scp, numpy as np
# from matplotlib import pyplot as plt
#
# import utils as uts
#
# PLOTS_PATH = 'plots/'
#
# if __name__ == '__main__':
#
#     # repeat N time:
#     # distribute R_all reads across l sequences uniformly and find the sequence with the minimal number of reads
#     # the plot is the histogram of the minimal values over N repeats
#     R_all = 4500
#     l = 50
#     N = 100000
#     T = 77
#     # Define a multinomial variable
#     X = scp.stats.multinomial(R_all,[1/l for _ in range(l)])
#     # generate N draws of X
#     x = X.rvs(N)
#     # Find N minimal values.
#     R_min = x.min(axis = 1)
#     #Plot distribution of minimal values
#     plt.figure(figsize=(10,6))
#     plt.subplots_adjust(bottom=0.2)
#     _ = plt.hist(R_min)#,bins = range(50,100))
#     plt.axvline(T,0,N/10,color='r')
#     # Calculate P(E)
#     P_E = np.sum(R_min < T) / N
#     # plt.legend(['T threshold', f'P(E)={P_E:.4f}', f'1-P(E)={1 - P_E:.4f}'])
#     plt.axvline(T, 0, N / 10, color='r', label='P(E) Threshold at T')
#     plt.axvline(R_all - T, 0, N / 10, color='g', label='1-P(E) Region Start')
#
#     plt.ylim(4000)
#     plt.title(f'$P(E)=P(min(R_j)<{T}) =$ {P_E}')
#     plt.xlabel('$min(R_j)$', fontsize=20)
#     plt.ylabel('Frequency', fontsize=20)
#     plt.yticks(fontsize=15)
#     plt.xticks(fontsize=15)
#     uts.make_dir(PLOTS_PATH)
#     plt.savefig(f"{PLOTS_PATH}/Multinomial bound R_all={R_all}, T={T}, l={l}, N={N}.svg", format='svg')

import scipy as scp, numpy as np
from matplotlib import pyplot as plt
import utils as uts

PLOTS_PATH = 'plots/'

if __name__ == '__main__':
    R_all = 5000
    l = 50
    N = 100000
    T = 77
    X = scp.stats.multinomial(R_all, [1/l for _ in range(l)])
    x = X.rvs(N)
    R_min = x.min(axis=1)

    plt.figure(figsize=(10,6))
    plt.subplots_adjust(bottom=0.2)
    counts, bins, patches = plt.hist(R_min, bins=range(min(R_min), max(R_min) + 1, 1), edgecolor='black')

    # Marking P(E) and 1-P(E)
    plt.axvline(T, color='r', label='P(E) Threshold at T')
    # plt.text(T+0.5, max(counts)/2, 'P(E)', color='red', verticalalignment='center')
    # plt.text(T-0.5, max(counts)/2, '1-P(E)', color='green', verticalalignment='center', horizontalalignment='right')

    P_E = np.sum(R_min < T) / N
    plt.title(f'Error Probability P(E) = {P_E}')
    plt.xlabel('Minimum Reads ($R_{\mathrm{min}}$)', fontsize=20)
    plt.ylabel('Frequency', fontsize=20)
    plt.legend()
    plt.yticks(fontsize=15)
    plt.xticks(fontsize=15)
    uts.make_dir(PLOTS_PATH)
    plt.savefig(f"{PLOTS_PATH}/Multinomial_bound_R_all={R_all}_T={T}_l={l}_N={N}.svg", format='svg')
    plt.show()
