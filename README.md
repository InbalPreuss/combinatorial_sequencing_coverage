# Combinatoric Coupon Collector

## Overview
The "Combinatoric Coupon Collector" project is designed to simulate and analyze a combinatorial version of the classic coupon collector's problem. This project is divided into three main parts:
1. **Reconstructing a Single Combinatorial Position**
2. **Reconstructing a Complete Combinatorial Sequence**
3. **Reconstructing a Complete Message**

Each part has both calculation and simulation modules, allowing for a comprehensive study of the problem from theoretical and empirical perspectives. Additionally, the project includes plotting scripts to visualize the outcomes of the simulations and calculations.

## Setup
To run this project, you will need Python and several libraries, including NumPy, Pandas, Matplotlib, and SciPy.

### Installation
1. Clone the repository:
   ```sh
   git clone [repository-url]

2. Install the required libraries using requirements.txt:
    ```sh
    pip install -r requirements.txt

## Usage
Navigate to the project directory and run the Python scripts for each part.

### Part 1: Single Position
* calculation/part1_reconstructing_a_single_combinatorial_position.py
* simulation/part1_reconstructing_a_single_combinatorial_position_simulation.py

Inputs:
Params in calculation and simulation:
* `n`: Total number of unique building blocks in each position.
* `t`: Required threshold on the number of observed occurrences.
* `R`: Acceptable error threshold.

Param only in simulation:
* `Q`: Number of iterations for the simulation to average the results and reduce variance.

Param only in calculation:
* `eps`: Error probability in the coupon collection process.

Output:
Probability distribution and visualization of reconstructing a single combinatorial position.

    python calculation/part1_reconstructing_a_single_combinatorial_position.py
    python simulation/part1_reconstructing_a_single_combinatorial_position_simulation.py

### Part 2: Complete Sequence
* calculation/part2_reconstructing_a_complete_combinatorial_sequence.py
* simulation/part2_reconstructing_a_complete_combinatorial_sequence_simulation.py

Inputs:
* Inherits `n`, `t`, `R` from Part 1.
* `m`: Total number of letters in the sequence.
* `b`: Number of letters required to be successfully reconstructed.

Param only in simulation:
* `Q`: Number of iterations for the simulation to average the results and reduce variance.

Param only in calculation:
* `eps`: Error probability in the coupon collection process.
* `method`: Statistical method used in calculations ('binomial' or 'normal').

Output:
Probability of successfully reconstructing the entire sequence.

    python calculation/part2_reconstructing_a_complete_combinatorial_sequence.py 
    python simulation/part2_reconstructing_a_complete_combinatorial_sequence_simulation.py

### Part 3: Complete Message
* calculation/part3_reconstructing_a_complete_message.py
* simulation/part3_reconstructing_a_complete_message_simulation.py

Inputs:
* Inherits `n`, `t`, `R`, `m`, `b`, from previous parts.
* `l`: Number of barcodes in the message.
* `a`: Number of barcodes required to be successfully decoded.

Param only in simulation:
* `Q`: Number of iterations for the simulation to average the results and reduce variance.

Param only in calculation:
* `eps`: Error probability in the coupon collection process.
* `delta`: Acceptable error threshold for the overall message decoding.
* `method`: Statistical method used in calculations ('binomial' or 'normal').

Output:
Decoding probability of the complete message.

    python calculation/part3_reconstructing_a_complete_message.py
    python simulation/part3_reconstructing_a_complete_message_simulation.py

### Running Main Simulations
Run main.py to execute simulations for all parts:

    python calculation/main.py

### Generating Plots
To generate plots, run plots.py:
   ```sh
   python plots.py
   ```


