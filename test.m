% M = cell2ndmat(h_seg);
% [~, n] = size(M);
% [D, T] = dist_seq(M(1:0.25*end, :).^2, 0.01, 5, 0.90, 5, 0.25);
% 
% R = T * ones(n, 1) / n;                                     % R: statistical results
% S = [];                                                     % S: resulting distance sequence
% for j = 1 : length(R)
%     if R(j) <= 0.25
%         S = [S; D(j, :)];
%     end
% end
% [b, a] = butter(2, 0.1);
% [len_S, ~] = size(S);
% sc = zeros(1, len_S);
% for i = 1 : len_S
%     sc(i) = test_autocorr(S(i, :), b, a);
% end
% sig = S(sc == min(sc), :);
% figure
% plot(sig)
% p = polyfit(1:length(sig), sig, 1);
% hold on
% plot(1:length(sig), p(1)*(1:length(sig))+p(2))

[m, n] = size(M);
h_ad = cell(1, n);
h_ad{1} = M(:, 1)';
for i = 2 : n
    h = M(:, i);
    hh = interp1(1:length(h), h, 1:0.01:length(h));
    step = -floor((i-1)*k*100);
    hh = [zeros(1, step), hh(1:end-step)];
    hhh = hh(1:100:end);
    h_ad{i} = hhh;
end
