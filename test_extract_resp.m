%% Find the most periodic range
[m, n] = size(M);
P = zeros(1, m);    
[b, a] = butter(2, 0.1);
for i = 1 : m
    P(i) = test_autocorr(M(i, :), b, a);
end
% figure
% stem(P)

% TODO: Extract the periodic range from P. Good luck :)
% ctr = find(P==max(P));
% range = ctr-20:ctr+20;
range = get_resp_range(P);

%% Extract the respiration curve
% figure
[D, T] = dist_seq(M(range, :), 0.1, 1, 0.3, 1, 0.25);
R = T * ones(n, 1) / n;                                     % R: statistical results
S = [];                                                     % S: resulting distance sequence
for j = 1 : length(R)
    if R(j) <= 0.25
        S = [S; D(j, :)];
    end
end
[len_S, ~] = size(S);
score = zeros(1, len_S);
for i = 1 : len_S
    score(i) = test_autocorr(S(i, :), b, a);
end
sig = S(score == max(score), :);
figure
plot(sig)
hold on
plot(filtfilt(b, a, sig))

%% Test

%% Backup
% for i = 1 : n 
%     l           = M(range, i);                              % One portion of h   
%     ll          = interp1(1:length(l), l, 1:0.1:length(l), 'spline');  % Interpolation
%     [pks, locs] = findpeaks(ll, 'MinPeakDistance', 30,...   % Find peaks & corresponding locations
%                             'MinPeakHeight', mean(ll));
%     if i == 1                                               % Initialization
%         D       = zeros(length(pks), n);                    % D: distance sequences
%         D(:, 1) = locs';
%         T       = zeros(length(pks), n);                    % T: testing sequences
%     else
%         for j = 1 : length(locs)                            % Iterate on each peak
%             flag = 0;                                       % flag: whether this peak is a new peak
%             for k = 1 : length(D(:, 1))                     % Iterate on each recorded peak
%                 if abs(locs(j) - D(k, i-1)) <= 10           % If the new peak is close to an existing peak, 
%                    D(k, i) = locs(j);                       % then they are the same peak
%                    flag    = 1;                             % This peak is clustered into one peak sequence
%                    break
%                 end
%             end
%             if flag == 0 && i <= 0.25*n                     % New peak is found & now is not too late
%                 D(length(D(:, 1))+1, i) = locs(j);          % Add a new peak sequence
%                 T(length(T(:, 1))+1, i) = 0;                % Corresponding testing sequence
%             end
%         end
%         for j = 1 : length(D(:, 1))                         % After all peaks are examined in this round,  
%             if D(j, i) == 0                                 % make up the sequences without new peak added
%                 D(j, i) = D(j, i-1);
%                 T(j, i) = 1;
%             end
%         end
%     end
% %     plot(ll) 
% %     hold on
% %     stem(locs, pks)
% %     hold off
% %     ylim([-0.4, 0.4])
% %     title(int2str(i))
% %     drawnow    
% end