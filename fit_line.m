function [k, b] = fit_line(h_seg, interp)
M = cell2ndmat(h_seg);
[~, n] = size(M);
[D, T] = dist_seq(M, 1:0.25*length(M(:, 1)), interp, 5, 0.90, 5, 0.25);

R = T * ones(n, 1) / n;                                     % R: statistical results
S = [];                                                     % S: resulting distance sequence
for j = 1 : length(R)
    if R(j) <= 0.25
        S = [S; D(j, :)];
    end
end
S = S*interp;
[len_S, ~] = size(S);
sc = zeros(1, len_S);
for i = 1 : len_S
    R = corrcoef(S(i, :), 1:length(S(i, :)));
    sc(i) = abs(R(1, 2));
end
sig = S(sc == max(sc), :);
p = polyfit(1:length(sig), sig, 1);
k = p(1);
b = p(2);
end