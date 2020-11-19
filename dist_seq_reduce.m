function [D, T] = dist_seq_reduce(D_last, M, range, int, peak_dis, peak_lvl, th, stp)
[~, n] = size(M);
D_last = D_last(:, 1:180);
D = D_last;
% figure
for i = 181 : n 
    l           = zeros(length(M(:, 1)));  
    l(range)    = M(range, i).^2;                              % One portion of h   
    ll          = interp1(1:length(l), l, 1:int:length(l), 'spline');  % Interpolation
    pk_h = prctile(ll, peak_lvl*100);
    [pks, locs] = findpeaks(ll, 'MinPeakDistance', peak_dis/int,...   % Find peaks & corresponding locations
                            'MinPeakHeight', pk_h);
%     plot(ll)
%     hold on
%     stem(locs, pks)
%     hold off
%     ylim([0 0.04])
%     drawnow
    if i == 1                                               % Initialization
        D       = zeros(length(pks), n);                    % D: distance sequences
        D(:, 1) = locs';
        T       = zeros(length(pks), n);                    % T: testing sequences
    else
        for j = 1 : length(locs)                            % Iterate on each peak
            flag = 0;                                       % flag: whether this peak is a new peak
            [mi, k] = min(abs(locs(j)-D(:,i-1)));           % Find the closest recorded peak
            if mi <= th/int                                 % Whether they belong to the same peak
                flag    = 1;
                D(k, i) = locs(j);
            end
            if flag == 0 && i <= stp*n                      % New peak is found & now is not too late
                D(length(D(:, 1))+1, i) = locs(j);          % Add a new peak sequence
                D(length(D(:, 1)), 1:i-1) = locs(j);
                T(length(T(:, 1))+1, i) = 0;                % Corresponding testing sequence
            end
        end
        for j = 1 : length(D(:, 1))                         % After all peaks are examined in this round,  
            if D(j, i) == 0                                 % make up the sequences without new peak added
                D(j, i) = D(j, i-1);
                T(j, i) = 1;
            end
        end
    end   
end
end