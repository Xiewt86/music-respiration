% =========================================================================
%
%           Coded by Wentao Xie - xiewt86@gmail.com
%
% =========================================================================
%% I. Constants & data preparation
clear
fs         = 48000;           % sampling rante of a smartphone
T_frame    = 0.1;             % sampling period of frames
T_fagg     = 0.4;             % compute CIR using T_fagg data 
L_frame    = fs*T_frame;      % length of a frame (CIR)
L_fagg     = fs*T_fagg;       % length of the aggregated data
T_agg      = 60;              % aggregate 60s of CIRs to compute resp rate
T_update   = 2;               % update the respiration rate every T_update
fs_e       = 1/T_frame;       % sampling rate of the CIR matrix
fs_g       = 20;              % sampling rate of the ground truth data
L_agg      = T_agg * fs;
L_agg_e    = T_agg * fs_e;    % # aggregated frames in estimation
L_update_e = T_update * fs_e; % # samples between each update in estimation
L_agg_g    = T_agg * fs_g;    % # aggregated data in ground truth data
L_update_g = T_update * fs_g; % # samples between each update of the ground truth data
music_code = 'orig0-';
echo_code  = 'XWT0_1';

% obtain coefficient k from .wav to .pcm
wav        = audioread('dataset/preamble.wav')'; %-1----clip of a .wav-----
pcm        = pcmread('dataset/preamble.pcm')'; %-1------clip of a .pcm-----
k          = max(pcm) / max(wav); %-1-------constant from .wav to .pcm-----
% 
org        = audioread(['dataset/', music_code, '.wav'])';
rcv        = pcmread(['dataset/raw-data/', echo_code, '.pcm'])' / k; %-1---
% rcv        = audioread(['dataset/raw-data/', echo_code, '.wav'])'; %-2-----

% synchronization: sb_music -> music
[f1, l1]   = xcorr(rcv(1:150*fs), org(1:150*fs));        % correlate rcv with org
[~, idx1]  = max(f1);                                    % match point
rcv        = rcv(l1(idx1) + 1 : l1(idx1) + length(org)); % correlation

%% II. Segmentation
% 1.low pass filter
[b, a]     = butter(20, 15000/(fs/2));
% 3.segmentation
offset     = 50;
len_cell   = ceil((length(org) - (L_fagg+offset) + 1) / (L_frame));
cnt_h = 1;
h_cell     = cell(1, len_cell);
cnt_sec    = ceil(length(org)/L_agg);
for i = 1 : cnt_sec
    if i < cnt_sec
        len_sec = L_agg;
    else
        len_sec = length(org) - (i-1)*L_agg;
    end
    len_h_tmp = ceil((len_sec - (L_fagg+offset) + 1) / (L_frame));
    h_cell_tmp = cell(1, len_h_tmp);
    for j = 1 : len_h_tmp
        x = org(ceil((i-1)*L_agg + 1+(j-1)*L_frame+offset) : ceil((i-1)*L_agg + 1+(j-1)*L_frame+offset+L_fagg-1));
        y = rcv((i-1)*L_agg + 1+(j-1)*L_frame : (i-1)*L_agg + 1+(j-1)*L_frame+L_fagg-1);
        h             = get_CIR(x, y, a, b);
        h_cell_tmp{j} = h(1 : 500);
    end
    M          = cell2ndmat(h_cell_tmp);
    F          = abs(range_fft(M)); 
    maxima     = max(F(:));
    [row, ~]   = find(F == maxima);
    if row - 50 >= 1
        low = row - 50;
    else
        low = 1;
    end
    if row + 50 <= 500
        high = row + 50;
    else
        high = 500;
    end
    for m = 1 : len_h_tmp
        h_cell_tmp{m} = h_cell_tmp{m}(low : high);
    end
    static_plot(h_cell_tmp);
    k = get_incline(h_cell_tmp);
    for j = 1 : len_h_tmp
        offset = offset + k;
        x      = org(ceil((i-1)*L_agg + 1+(j-1)*L_frame+offset) : ceil((i-1)*L_agg + 1+(j-1)*L_frame+offset+L_fagg-1));
        y      = rcv((i-1)*L_agg + 1+(j-1)*L_frame : (i-1)*L_agg + 1+(j-1)*L_frame+L_fagg-1);
        h             = get_CIR(x, y, a, b);
        h_cell{cnt_h} = h(1 : 500);
        cnt_h         = cnt_h+1;
    end
end

for i = 1 : length(h_cell)
    if isempty(h_cell{i})
        h_cell = h_cell(1:i-1);
        break
    end
end
        
% 4.eliminate the odd h (hugely different from others)
for i = 2 : length(h_cell)-1
    if sum(abs(h_cell{i} - h_cell{i-1})) + sum(abs(h_cell{i} - h_cell{i+1})) > 20
        h_cell{i} = (h_cell{i-1} + h_cell{i+1})/2;
    end
end


%% III. Analysis             
ranges        = struct();        % the analysis ranges of the range-frequency matrix
ranges.rBins  = 10 : 150;         % ... range bin
ranges.fBins  = 10 : 30;         % ... frequency bin

rate_len      = ceil((length(h_cell)-L_agg_e+1)/L_update_e);   % total length
rates_e       = zeros(1, rate_len);                            % estimated respiration rate
rates_g       = zeros(1, rate_len);                            % ground truth respiration rate

preamble      = 60.1;            % length of the preamble (in sec)
delay         = 20;              % offset produced manually
data_struct   = load(['dataset/ground-truth/', echo_code, '.mat']);
gt_raw        = data_struct.data(preamble*fs_g+delay : end);

for i = 1 : rate_len
    h_seg      = h_cell(L_update_e*(i-1)+1 : L_update_e*(i-1) + L_agg_e);  % aggregate T_agg data from CIRs
    gt_seg     = gt_raw(L_update_g*(i-1)+1 : L_update_g*(i-1)+L_agg_g);    % aggregate T_agg data from ground truth
    if rates_g(1) == 0
        [rate, bins] = estimate_frequency(h_seg, -1, -1, ranges);
    else
        [rate, bins] = estimate_frequency(h_seg, rates_e(i-1), bins, ranges);
    end
    rates_e(i) = rate;
    rates_g(i) = ground_frequency(gt_seg, ranges);
end

figure
plot(rates_e)
hold on
plot(rates_g)
static_plot(h_cell)

%% -------------------function 1 get_frequency ----------------------------
function [f, bins] = estimate_frequency(h_cell, f_last, b_last, setting)
M          = cell2ndmat(h_cell);                   % CIR matrix
F          = abs(range_fft(M));                    % frequency matrix
F_valid    = F(setting.rBins, setting.fBins);      % target frequency matrix
maxima     = max(F_valid(:));
[row, ~]   = find(F_valid == maxima);
f_vector   = F_valid(row, :);
[f, f_vec] = vector_estimate(f_vector, maxima, 0.7, setting);
bins       = f_vec;
if f_last == -1
    return
end
if abs(f-f_last) > 3
    f_vec = b_last;
    w_vec = f_vector(b_last);
    f     = f_vec*w_vec'/sum(w_vec);
    bins  = f_vec;
end
if abs(f-f_last) > 3
    f    = f_last;
    bins = f_vec;
end
end

%% --------------------- function 2 get_CIR -------------------------------
function h = get_CIR(x, y, a, b)
X     = fft(x);
k0    = X==0;
X(k0) = 0.0001;
Y     = fft(y);
H     = Y ./ X;
h_raw = real(ifft(H));
h_tmp = mapminmax(filter(b, a, h_raw));
h     = h_tmp - mean(h_tmp);
end

%% -------------------- function 3 cell2ndmat -----------------------------
function M = cell2ndmat(hh)
n = length(hh);
m = length(hh{1});
M = zeros(m, n);
for i = 1 : n
   M(:, i) = hh{i}';
end
end

%% --------------------- function 4 range_fft -----------------------------
function F = range_fft(M)
[m, n] = size(M);
F      = zeros(m, n);
for i = 1 : m
   F(i, :) = (fft(M(i, :))); 
end
end

%% ------------------ function 5 ground_frequency -------------------------
function F = ground_frequency (gt, setting)
GT     = abs(fft(gt));
GT     = GT(setting.fBins);
maxima = max(GT);
[F, ~] = vector_estimate(GT, maxima, 0.7, setting);
end

%% ------------------- function 6 vector_estimate -------------------------
function [f, f_vec] = vector_estimate(F_vector, maxima, Th, setting)
f_vec  = [];
w_vec  = [];
for j = 1 : length(F_vector)
    if F_vector(j) > Th*maxima
        w_vec = [w_vec, F_vector(j)];
        f_vec = [f_vec, j+setting.fBins(1)-2];
    end
end
f = f_vec*w_vec'/sum(w_vec);
end
%%
function k = get_incline(h_cell)
h_model = h_cell{1};
idx     = zeros(1, length(h_cell));
for i = 1 : length(h_cell)
    [f, l] = xcorr(h_model, h_cell{i});
    [~, j] = max(f);
    idx(i) = l(j);
end
[~, locs] = findpeaks(-idx, 'MinPeakHeight', -mean(idx));
idx2       = idx(locs);
[p, ~]    = polyfit(locs, idx2, 1);
k         = p(1);
end
