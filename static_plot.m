function M = static_plot(C)
M = [];
f1 = figure;
for i = 1 : length(C)
    M = [M, C{i}'];
end
imagesc(M)
colormap(f1, gray)
end
