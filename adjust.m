h_model = h_cell{1};
idx = zeros(1, length(h_cell));
for i = 1 : length(h_cell)
    [f, l] = xcorr(h_model, h_cell{i});
    [~, j] = max(f);
    idx(i) = l(j);
end
figure
plot(idx)
[~, locs] = findpeaks(-idx, 'MinPeakHeight', -mean(idx));
% idx = idx(locs);
[p, s] = polyfit(locs, idx(locs), 1);
hold on
plot(locs, idx(locs))
hold on
plot(1:length(idx), p(1)*(1:length(idx)) + p(2))