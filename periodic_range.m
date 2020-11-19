function [range, ctr] = periodic_range(M, last_ctr) 
if last_ctr == -1
    [m, ~] = size(M);
    itr = 1 : m;
else
    m = last_ctr;
    itr = m-30 : m+30;
end
[b, a] = butter(2, 0.1);
for i = itr
    P(i) = test_autocorr(M(i, :), b, a);
end

% TODO: Extract the periodic range from P. Good luck :)
ctr = find(P==max(P));
range = ctr-20:ctr+20;
% range = get_resp_range(P);
end
