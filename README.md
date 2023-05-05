# Closed-form Weight Enumeration of Polar Codes
This MATLAB script enumerates the minimum weight and 1.5 times minimum weight codewords of decreasing monomial codes including polar codes and Read-Muller in no time and returns the followings: 

$w = [ w_{min}, 1.5w_{min} ]$

$A_w = [ A_{w_{min}}, A_{1.5w_{min}} ]$

These two weights are the smallest weights of polar codes. These two weights are probably the most important weights as they are dominant in computing the union bound (the upper bound for block error rate). 

The algorithm is super fast because the computation is based on closed form expressions.

## Inputs:
- I: Indices of $K$-most relaible bit-channels, allocated for information bits
- n: $\log_2N$ where $N$ is the code length

## Outputs: 
- r: the magimum degree of monomials. minimum distance = $2^{n-r}$. 
- w: a vector indicating the minimum weight $w_{min}$ and $1.5w_{min}$ of the code.
- A_w: a vector containing the multiplicities of codewords with weights w_min and $1.5w_{min}$, correspoding to output $w$.
- A: a matrix showing the breakdown of cosets/subgroups generating 1.5w_min-weight codewords. For details, see Table I in https://arxiv.org/abs/2305.02921

---
This script showcases the results in the paper below:

M. Rowshan, V.F. Dragoi, and J. Yuan, “On the Closed-form Weight Enumeration of Polar Codes: 1.5d-weight Codewords”, preprint, 2023. arXiv:2305.02921

Please report any bugs to mrowshan at ieee dot org
