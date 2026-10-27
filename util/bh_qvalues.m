function q = bh_qvalues(p)
%BH_QVALUES Benjamini-Hochberg false discovery rate q values.
%
% q = bh_qvalues(p) returns, for each p value, the smallest false
% discovery rate at which it would be declared significant under the
% Benjamini and Hochberg (1995) procedure. p must not contain NaN.
% Equivalent to mafdr(p,'BHFDR',true) without the Bioinformatics Toolbox.
%
% Alex Kinsella, October 2026

sz = size(p);
p = p(:);
m = numel(p);

[p_sorted, ord] = sort(p);
q_sorted = p_sorted .* (m ./ (1:m)');
q_sorted = min(cummin(q_sorted, 'reverse'), 1);

q = nan(m,1);
q(ord) = q_sorted;
q = reshape(q, sz);
end
