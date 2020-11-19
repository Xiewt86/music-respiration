function [D, T] = dist_seq(M, range, int, peak_dis, peak_lvl, th, stp)
[~, n] = size(M);
% figure
for i = 1 : n 
    l           = M(range, i);                              % One portion of h   
    ll          = interp1(1:length(l), l, 1:int:length(l), 'spline');  % Interpolation
    pk_h = prctile(ll, peak_lvl*100);
    [pks, locs] = findpeaks(ll, 'MinPeakDistance', peak_dis/int,...   % Find peaks & corresponding locations
                            'MinPeakHeight', pk_h);
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
%     plot(ll)
%     hold on
%     stem(D(:, i), ll(D(:, i)))
%     hold off
%     ylim([-0.4 0.4])
%     drawnow
%         pause 
end
D = D+range(1);
end