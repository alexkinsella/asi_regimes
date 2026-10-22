function [p_cutoff,cnt] = BH(pval,alpha)

% Benjamini-Hochberg Procedure
% Control the false discovery rate within correlations
% Input: p-values and control level alpha
% Output: Cutoff p-value and its index in the sorted values
m = numel(pval);

[p_srt,srt_ind] = sort(pval);

p_good = nan(size(p_srt));
for kk = 1:m
    if p_srt(kk)>kk*alpha/m
        p_good(kk) = 0;
    else
        p_good(kk) = 1;
    end
end
cnt = find(p_good,1,'last');
p_cutoff = p_srt(cnt);