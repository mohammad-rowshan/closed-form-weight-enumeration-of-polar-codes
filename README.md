# Closed-form Weight Enumeration of Polar Codes
The MATLAB script enumerates the minimum weight and 1.5 times minimum weight codewords of polar codes. These two weights are the smallest weights of polar codes. The algorithm is super fast because the computation is based on closed form expressions.

#Inputs:
I: Indices of K-most relaible bit-channels, allocated for information bits

n: log2(N) where N is the code length

#Outputs: 
r: the magimum degree of monoials. minimum distance = 2^(n-r). 
w: a vector indicating the minimum weight (w_min) and 1.5w_min.

A_w: a vector containing the multiplicities of codewords with weights w_min and 1.5w_min, correspoding to output w.

A: a matrix showing the breakdown of cosets/subgroups generating 1.5w_min-weight codewords. For details, see Table I in https://arxiv.org/abs/2305.02921

This script showcases the results in the paper below:

M. Rowshan, Vlad-Florin Dragoi, and Jinhong Yuan, “On the Closed-form Weight Enumeration of Polar Codes: 1.5d-weight Codewords”, preprint, 2023. arXiv:2305.02921

Please report any bugs to mrowshan at ieee dot org
