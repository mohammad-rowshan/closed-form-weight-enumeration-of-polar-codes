% Weight Enumeration of Polar Codes: First two smallest weights %%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (c) 2023, Mohammad Rowshan, Vlad-Florin Dragoi, and Jinhong Yuan
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without modification, 
% are permitted provided that:
% the source code retains the above copyright notice, and the redistribtuion condition.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = 7
N = 2^n;            % Code length
R = 0.5;            % Code rate
K = N * R;          % code dimention: the number of information bits
design_snr_db = 0;  % for the optimized code

I = construct_dega(design_snr_db, N, K)'; % the indices of K information bits

[r,w,A_w,A] = enum_onehalf_w(sort(I),n) 

% returns the minimum distance of the code and the number of codewords with minumum distance (error coefficient)
function [r,w, A_w, A] = enum_onehalf_w(I,n)
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
                A = [A; [i,j,h-1,foh-1,goh-1,alpha,A_fg]];
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



%%%%%%%%%%%%% Code Constructin %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% retruns Mean-LLRs (a measure for sub-channels' reliability) obtained from Density Evolution by Gaussian Approximation (DEGA) by Chung et. al. and Trifonov, Algorithm by Vangala et. al.
function I =  construct_dega(design_snr_db, N, K)
    mllr = zeros(N,1);
    sigma_sq = 1/(2*K/N*power(10,design_snr_db/10));
    mllr(1) = 2/sigma_sq;
    for level = 1: log2(N)
        B = 2^level;
        for j = 1 :B / 2
            T = mllr(j);
            mllr(j) = calc_phi_inv(T);
            mllr(B / 2 + j) = 2 * T;
        end
    end
    
    mask = zeros(N,3);
    for i = 0:N-1
        nat(i+1) = bitreversed(i,uint8(log2(N)));
    end
    %nat = bitrevorder(0:N-1);
    for i = 1:N
        mask(i,:) = [nat(i), mllr(i), 1];
    end
    % sort sub-channels by mllr
    mask = sortrows(mask,2); %direction: ascend (default)
    % set info bits to 1 for sub-channels with K largest mllr values
    for i = 1:N-K
        mask(i,3) = 0;
    end
    % sort channels with respect to index (in bitreversal order; line 42
    mask = sortrows(mask,1); %direction: ascend (default)
    I = find(mask(:,3)==1)-1;
end

function dec = bitreversed(i,n) % You can instead use bitrevorder() in "Singal Processing" toolbox.
    dec = bin2dec(fliplr(dec2bin(i,n)));
end

% returns Phi inverse based on piece-wise linear approximation, by Trifonov
function phi_inv = calc_phi_inv(x)
    if (x>12)
        phi_inv = 0.9861 * x - 2.3152;
    elseif (x<=12 && x>3.5)
        phi_inv = x*(0.009005 * x + 0.7694) - 0.9507;
    elseif (x<=3.5 && x>1)
        phi_inv = x*(0.062883*x + 0.3678)- 0.1627;
    else
        phi_inv = x*(0.2202*x + 0.06448);
    end
end
