from math import ceil

from part1_reconstructing_a_single_combinatorial_position import ReconstructingSingleCombinatorialPosition
from part2_reconstructing_a_complete_combinatorial_sequence import ReconstructingCombinatorialSequence
from part3_reconstructing_a_complete_message import MessageDecoder


# Press the green button in the gutter to run the script.
if __name__ == '__main__':

    ##########
    # PART 1 #
    ##########

    n = 7  # Total number of unique building blocks in each position
    t = 3  # Required threshold on the number of observed occurrences
    eps = 0.01
    R = 50  # Acceptable error threshold

    # Run
    part1_reconstructing_single_combinatorial_position = ReconstructingSingleCombinatorialPosition(n, t, eps, R)
    rr, p, r = part1_reconstructing_single_combinatorial_position.calc_results()

    ##########
    # PART 2 #
    ##########

    # Part 1 params
    n = 5  # Total number of unique building blocks in each position
    t = 3  # Required threshold on the number of observed occurrences
    eps = 0.01
    R = 40  # Acceptable error threshold

    # Part 2 params
    m = 10  # Total number of letters in the sequence
    b = 9  # Number of letters required to be successfully reconstructed
    method = "binomial"

    # Run
    reconstructor = ReconstructingCombinatorialSequence(n=n, t=t, eps=eps, R=R, m=m, b=b, method=method)
    P_succ_single, pi_R = reconstructor.calculate_probability()

    ##########
    # PART 3 #
    ##########

    # Part 1 params
    n = 5  # Total number of unique building blocks in each position
    t = 3  # Required threshold on the number of observed occurrences
    eps = 0.01
    R = 40  # Acceptable error threshold

    # Part 2 params
    m = 120  # Total number of letters in the sequence
    b = 100  # Number of letters required to be successfully reconstructed
    method = "binomial"

    # Part 3 params
    l = 10  # Number of barcodes in the message
    a = ceil(l * 0.8)  # Number of barcodes required to be successfully decoded
    # eps = 0.01
    delta = 0.1  # Acceptable error threshold
    P_E_method = 'simulation'
    # Run
    decoder = MessageDecoder(n=n, t=t, eps=eps, R=R, m=m, b=b, method=method, l=l, a=a, delta=delta, P_E_method=P_E_method)
    R_all = decoder.run_decoder()

