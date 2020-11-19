function M = cell2ndmat(hh)
n = length(hh);
m = length(hh{1});
M = zeros(m, n);
for i = 1 : n
   M(:, i) = hh{i}';
end
end