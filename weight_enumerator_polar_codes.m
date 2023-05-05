% Weight Enumeration of Polar Codes: First two smallest weights %%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (c) 2023, Mohammad Rowshan, Vlad-Florin Dragoi, and Jinhong Yuan
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without modification, 
% are permitted provided that:
% the source code retains the above copyright notice, and the redistribtuion condition.
%
% weight_enum function Inputs/Outputs%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
% I: Indices of K-most relaible bit-channels, allocated for information bits
% n: log2(N) where N is the code length
% Outputs:
% r: the magimum degree of monoials. minimum distance = 2^(n-r). 
% w: a vector indicating the minimum weight (w_min) and 1.5*w_min.
% A_w: a vector containing the multiplicities of codewords with weights
% w_min and 1.5*w_min, correspoding to output w.
% A: a matrix showing the breakdown of cosets/subgroups generating
% 1.5w_min-weight codewords. For details, see Table I in https://arxiv.org/abs/2305.02921
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = 7
N = 2^n;            % Code length
R = 0.5;            % Code rate
K = N * R;          % code dimention: the number of information bits

design_snr_db = 0;  % for the optimized code
I = rate_profile(design_snr_db, N, K)'; % the indices of K information bits

[r,w,A_w,A] = weight_enum(sort(I),n)

% returns the minimum distance of the code and the number of codewords with minumum distance (error coefficient)
function [r,w, A_w, A] = weight_enum(I,n)
    r = max(sum(~(dec2bin(I)-'0'),2));
    Ir = I(find(sum(~(dec2bin(I)-'0'),2)==r)); % The ascending order in I is needed for conditions in computing alpha
    Ir_sub = Ir; w=zeros(1,2); A_w=zeros(1,2);  w(1)=2^(n-r); w(2)=1.5*w(1); A = [];
    for i = Ir
        f = find(~(reverse(dec2bin(i,n))-'0'));
        A_w(1) = A_w(1) + 2^(r+lambda(f,f));
        Ir_sub(Ir_sub==i) = [];
        for j = Ir_sub
            g = find(~(reverse(dec2bin(j,n))-'0'));
            h = intersect(f,g);
            if length(h)==r-2
                foh = setdiff(f,h); goh = setdiff(g,h);
                alpha = (foh(2)>goh(2) & goh(2)>foh(1)) + (goh(1)>foh(1));
                A_fg = 2^(r+2 + lambda(h,h) + lambda(f,foh) + lambda(g,goh) - alpha);
                A_w(2) = A_w(2) + A_fg; 
                A = [A; [i,j,h-1,foh-1,goh-1,alpha,A_fg]]; % Provides
            end
        end
    end
end
function orbit =  lambda(f,g)
    orbit = 0;
    for i = g
        orbit = orbit + length(setdiff([1:i-1],f));
    end
end
