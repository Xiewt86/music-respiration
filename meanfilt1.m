function s = meanfilt1(sig, win)
s = zeros(1, length(win)-floor(win/2)+1);
for i = floor(win)+1 : length(sig)-floor(win/2)+1
    s(i) = mean(sig(i-floor(win/2) : i+win-floor(win/2)-1));
end
end

