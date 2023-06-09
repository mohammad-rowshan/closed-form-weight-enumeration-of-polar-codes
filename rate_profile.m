% retruns the indices of K-most relaible bit-channels, set I.
% The measure for bit-channels' reliability is Mean-LLRs obtained from Density Evolution 
% by Gaussian Approximation (DEGA) by Chung et. al. and Trifonov, Algorithm by Vangala et. al.
function I =  rate_profile(design_snr_db, N, K)
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