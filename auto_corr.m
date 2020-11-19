function f = auto_corr(x)
f = [];
for i = 1 : length(x)
    sum = 0;
    for j = i : length(x)
        sum = sum + x(j)*conj(x(j-i+1));
    end
    f = [f, sum];
end
