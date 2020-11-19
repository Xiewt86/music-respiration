function l = line_plot(C, raw)
M = [];
for i = 1 : length(C)
    M = [M, C{i}'];
end
% figure
% plot(M(raw, :))
l = M(raw, :);
end
