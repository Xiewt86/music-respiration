function [f, l] = get_rate(curve, fs, Th)
C = abs(fft(curve, 10000));
C = C(1:floor(6667/fs));
[pks, locs] = findpeaks(C);
pks = pks(2:end);
locs = locs(2:end);
M = max(pks);
l = [];
for i = 1 : length(locs)
    if pks(i) >= Th*M
        l = [l, locs(i)];
    end
end
f = (l*C(l)'/sum(C(l))-1)*fs/10000*60;
end