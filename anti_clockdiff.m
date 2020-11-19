function h_ad = anti_clockdiff(h_seg, k, interp)
n = length(h_seg);
h_ad = cell(1, n);
h_ad{1} = h_seg{1};
for i = 2 : n
    h = h_seg{i};
    hh = interp1(1:length(h), h, 1:interp:length(h));
    step = -floor((i-1)*k/interp);
    hh = [zeros(1, step), hh(1:end-step)];
    hs = hh(1:1/interp:end);
    h_ad{i} = hs;
end
end