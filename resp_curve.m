function [curve, f, S] = resp_curve(M, range, b, a)
f = 0;
[~, n] = size(M);
[D, T] = dist_seq(M, range, 0.1, 1, 0.5, 2, 0.25);       % D: Distance sequences
R = T * ones(n, 1) / n;                                     % R: statistical results
S = [];                                                     % S: resulting distance sequence
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
if isempty(score) || max(score) < 100
    f = 1;
    score = zeros(1, length(M(:, 1)));
    for i = range
        score(i) = test_autocorr(M(i, :), b, a);
    end
    sig = -M(score == max(score), :);
else
    sig = S(score == max(score), :);
end
curve = filtfilt(b, a, sig);
end
