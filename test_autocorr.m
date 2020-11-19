function delta = test_autocorr(sig, b, a)
f = auto_corr(zscore(sig));
ff = filtfilt(b, a , f);
[~, locs1] = findpeaks(ff);
[~, locs2] = findpeaks(-ff);

% sum = 0;
% i_p = 1;
% i_v = 1;
% while i_p <= length(locs1)
%     if locs1(i_p) <= locs2(i_v)
%         i_p = i_p + 1;
%     elseif i_v >= length(locs2)
%         sum = sum + ff(locs1(i_p))-ff(locs2(i_v));
%         break
%     elseif locs2(i_v+1) <= locs1(i_p)
%         i_v = i_v + 1;
%     else 
%         sum = sum + ff(locs1(i_p))-ff(locs2(i_v));
%         i_v = i_v + 1;
%     end
% end
% delta = sum;

if isempty(locs1) || isempty(locs2)
    delta = 0;
else
    delta = ff(locs1(1))-ff(locs2(1));
end
end