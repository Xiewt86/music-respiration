clear
fs         = 48000;       % sampling rante of a smartphone
T_frame    = 0.1;         % sampling period of frames  0.1 by default
T_fagg     = 0.4;         % compute CIR using T_fagg data    0.4 by default
L_frame    = fs*T_frame;  % length of a frame (CIR)
L_fagg     = fs*T_fagg;   % length of the aggregated data
music_code = 'orig1';

echo_code  = 'XWT1';
clock_diff = 0;
step_show  = 0;

% obtain coefficient k from .wav to .pcm
wav        = audioread('dataset/preamble.wav')'; %-1----clip of a .wav-----
pcm        = pcmread('dataset/preamble.pcm')'; %-1------clip of a .pcm-----
k          = max(pcm) / max(wav); %-1-------constant from .wav to .pcm-----

org        = audioread(['dataset2/', music_code, '.wav'])';
rcv        = pcmread(['dataset2/raw-data/', echo_code, '.pcm'])' / k; %-1---
% rcv        = audioread(['dataset/raw-data/', echo_code, '.wav'])'; %-2-----

% synchronization: sb_music -> music
[f1, l1]   = xcorr(rcv(1:150*fs), org(1:150*fs));        % correlate rcv with org
[~, idx1]  = max(f1);                                    % match point
rcv        = rcv(l1(idx1) + 1 : l1(idx1) + length(org)); % correlation

%% II. Segmentation
% 1.low pass filter
[b, a]     = butter(10, 15000/(fs/2));
% 2.find the index of the maxima in CIR
h_tmp      = zeros(1, L_fagg);
for i = 1 : L_frame : 1+50*L_frame
    x     = org(i+50 : i+L_fagg-1+50);
    y     = rcv(i : i+L_fagg-1);
    h_tmp = h_tmp + get_CIR(x, y, a, b);
end
maxima     = max(h_tmp);
dmax       = find(h_tmp == maxima);           % a CIR should start at 'dmax'
% dmax = 30;  %---------------------------------------------------------------
% 3.segmentation
offset     = 50;
cell_len   = ceil((length(org) - (L_fagg+offset) + 1) / (L_frame));
h_cell     = cell(1, cell_len);
for i = 1 : cell_len
    x         = org(1+(i-1)*L_frame+offset : 1+(i-1)*L_frame+offset+L_fagg-1);
    y         = rcv(1+(i-1)*L_frame : 1+(i-1)*L_frame+L_fagg-1);
    h         = get_CIR(x, y, a, b);
    if sum(h(dmax : dmax+599).^2)>=20 && i >= 100
        h_cell{i} = h_cell{i-1};
    else
        h_cell{i} = h(dmax : dmax+599);
    end
end

%% III. Analysis             
T_agg         = 20;              % aggregate 60s of CIRs to compute resp rate
T_update      = 2;               % update the respiration rate every T_update
fs_e          = 1/T_frame;       % sampling rate of the CIR matrix
fs_g          = 20;              % sampling rate of the ground truth data
L_agg_e       = T_agg * fs_e;    % # aggregated frames in estimation
L_update_e    = T_update * fs_e; % # samples between each update in estimation
L_agg_g       = T_agg * fs_g;    % # aggregated data in ground truth data
L_update_g    = T_update * fs_g; % # samples between each update of the ground truth data

% rate_len      = ceil((length(h_cell)-L_agg_e+1)/L_update_e);   % total length
% rates_e       = zeros(1, rate_len);                            % estimated respiration rate
% rates_g       = zeros(1, rate_len);                            % ground truth respiration rate

[B, A] = butter(2, 0.1); 
last_ctr = -1;
[curve, f, last_ctr] = get_resp(h_cell, B, A, last_ctr);
e = zscore(curve);

preamble      = 17.9;            % length of the preamble (in sec)
delay         = 1*fs_g;              % offset produced manually
data_struct   = load(['dataset2/ground-truth/', echo_code, '.mat']);
gt_raw        = data_struct.data(preamble*fs_g+delay : end)';
gt            = -gt_raw(1:2:end);
g = zscore(gt(1:length(e)));

figure
plot(e)
hold on
plot(g)
hold off

rate_len = ceil((length(h_cell)-L_agg_e+1)/L_update_e);   % total length
rates_e = []; rates_g = [];
rec = [];
for i = 1 : rate_len
    disp([int2str(i), ' out of ', int2str(rate_len)])
    e_seg = e(L_update_e*(i-1)+1 : L_update_e*(i-1) + L_agg_e);
    g_seg = g(L_update_e*(i-1)+1 : L_update_e*(i-1) + L_agg_e);
    rates_g = [rates_g, get_rate(g_seg, fs_e, 0.9)];
    rates_e = [rates_e, get_rate(e_seg, fs_e, 0.9)];
%     e_r = get_rate(e_seg, fs_e, 0.9);
%     if abs(e_r-rates_e) > 5
%         rates_e = [rates_e, rates_e(end)];
%     else
%         rates_e = [rates_e, e_r];
%     end
end

figure
plot(rates_e)
hold on
plot(rates_g)
hold off

% function f = get_rate(curve, fs, Th)
% C = abs(fft(curve, 10000));
% C = C(1:floor(6667/fs));
% [pks, locs] = findpeaks(C);
% pks = pks(2:end);
% locs = locs(2:end);
% M = max(pks);
% l = [];
% for i = 1 : length(locs)
%     if pks(i) >= Th*M
%         l = [l, locs(i)];
%     end
% end
% f = (l*C(l)'/sum(C(l))-1)*fs/10000*60;
% end

function [curve, f, ctr] = get_resp(h_seg, B, A, last_ctr)
M = cell2ndmat(h_seg);
[range, ctr] = periodic_range(M, last_ctr);
[curve, f] = resp_curve(M, range, B, A);
'Pause Here';
end
