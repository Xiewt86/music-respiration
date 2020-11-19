function [curve, f, S] = resp_curve_reduce(M, range, b, a, D_last)
[~, n] = size(M);
[D, T] = dist_seq_reduce(D_last, M, range, 0.1, 1, 0.5, 1, 0.25);       % D: Distance sequences
R = T * ones(n, 1) / n;                                                 % R: statistical results
S = [];                                                                 % S: resulting distance sequence
for j = 1 : length(R)
    if R(j) <= 0.25
        S = [S; D(j, :)];
    end
end

[len_S, ~] = size(S);
score = zeros(1, len_S);
for i = 1 : len_S
    score(i) = test_autocorr(S(i, :), b, a);
end
if max(score) >= 300
    sig = S(score == max(score), :);
    curve = filtfilt(b, a, sig);
    f = 0;
else
    f = 1;
end
