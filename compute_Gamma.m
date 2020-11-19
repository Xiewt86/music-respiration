% =========================================================================
%
%           Coded by Wentao Xie - xiewt86@gmail.com
%
% =========================================================================
%% 
clear
fs         = 48000;       % sampling rante of a smartphone
T_frame    = 0.1;         % sampling period of frames  0.1 by default
T_fagg     = 0.4;         % compute CIR using T_fagg data    0.4 by default
L_frame    = fs*T_frame;  % length of a frame (CIR)
L_fagg     = fs*T_fagg;   % length of the aggregated data
music_code = 'orig0-';
echo_code  = 'XWT0_1';

% obtain coefficient k from .wav to .pcm
wav        = audioread('dataset/preamble.wav')'; %-1----clip of a .wav-----
pcm        = pcmread('dataset/preamble.pcm')'; %-1------clip of a .pcm-----
k          = max(pcm) / max(wav); %-1-------constant from .wav to .pcm-----

org        = audioread(['dataset/', music_code, '.wav'])';
rcv        = pcmread(['dataset/raw-data/', echo_code, '.pcm'])' / k; %-1---
% rcv        = audioread(['dataset/raw-data/', echo_code, '.wav'])'; %-2-----

% synchronization: sb_music -> music
[f1, l1]   = xcorr(rcv(1:150*fs), org(1:150*fs));        % correlate rcv with org
[~, idx1]  = max(f1);                                    % match point
rcv        = rcv(l1(idx1) + 1 : l1(idx1) + length(org)); % correlation

% 1.low pass filter
[b, a]     = butter(20, 15000/(fs/2));
% 2.find the index of the maxima in CIR
h_tmp      = zeros(1, L_fagg);
for i = 1 : L_frame : 1+50*L_frame
    x     = org(i+50 : i+L_fagg-1+50);
    y     = rcv(i : i+L_fagg-1);
    h_tmp = h_tmp + get_CIR(x, y, a, b);
end
maxima     = max(h_tmp);
dmax       = find(h_tmp == maxima);           % a CIR should start at 'dmax'
% dmax = 1;  %---------------------------------------------------------------
% 3.segmentation
offset     = 50;
cell_len   = ceil((length(org) - (L_fagg+offset) + 1) / (L_frame));
h_cell     = cell(1, cell_len);
record = zeros(1, cell_len);
for i = 1 : cell_len
    x         = org(1+(i-1)*L_frame+offset : 1+(i-1)*L_frame+offset+L_fagg-1);
    y         = rcv(1+(i-1)*L_frame : 1+(i-1)*L_frame+L_fagg-1);
    h         = get_CIR(x, y, a, b);
    h_cell{i} = h;
end

Gamma_cell4 = {};
N = length(h_cell{1});
tal = 1129;
for cnt = 1 : 10 : length(h_cell)-1
    disp(['h_cell no. ', int2str(cnt), ' total ', int2str(length(h_cell)-1)]);
    x = h_cell{cnt};
    xp = h_cell{cnt+1};
    Gamma = zeros(1, N);
    for k = 0 : N-1
        common = 0;
        for n = 0 : N-tal-1
            common = common + x(n+1)*exp(-1j*2*pi*k*n/N);
        end
        lower = common;
        for n = N-tal : N-1
            lower = lower + x(n+1)*exp(-1j*2*pi*k*n/N);
        end
        upper = common;
        for n = N-tal : N-1
            upper = upper + xp(n+1)*exp(-1j*2*pi*k*n/N);
        end
        Gamma(k+1) = upper / lower;
    end
    Gamma_cell4 = [Gamma_cell4, Gamma];
end

%%
fs         = 48000;       % sampling rante of a smartphone
T_frame    = 0.1;         % sampling period of frames  0.1 by default
T_fagg     = 0.3;         % compute CIR using T_fagg data    0.4 by default
L_frame    = fs*T_frame;  % length of a frame (CIR)
L_fagg     = fs*T_fagg;   % length of the aggregated data
music_code = 'orig0-';
echo_code  = 'XWT0_1';

% obtain coefficient k from .wav to .pcm
wav        = audioread('dataset/preamble.wav')'; %-1----clip of a .wav-----
pcm        = pcmread('dataset/preamble.pcm')'; %-1------clip of a .pcm-----

org        = audioread(['dataset/', music_code, '.wav'])';
rcv        = pcmread(['dataset/raw-data/', echo_code, '.pcm'])' / k; %-1---
% rcv        = audioread(['dataset/raw-data/', echo_code, '.wav'])'; %-2-----

% synchronization: sb_music -> music
[f1, l1]   = xcorr(rcv(1:150*fs), org(1:150*fs));        % correlate rcv with org
[~, idx1]  = max(f1);                                    % match point
rcv        = rcv(l1(idx1) + 1 : l1(idx1) + length(org)); % correlation

% 1.low pass filter
[b, a]     = butter(20, 15000/(fs/2));
% 2.find the index of the maxima in CIR
h_tmp      = zeros(1, L_fagg);
for i = 1 : L_frame : 1+50*L_frame
    x     = org(i+50 : i+L_fagg-1+50);
    y     = rcv(i : i+L_fagg-1);
    h_tmp = h_tmp + get_CIR(x, y, a, b);
end
maxima     = max(h_tmp);
dmax       = find(h_tmp == maxima);           % a CIR should start at 'dmax'
% dmax = 1;  %---------------------------------------------------------------
% 3.segmentation
offset     = 50;
cell_len   = ceil((length(org) - (L_fagg+offset) + 1) / (L_frame));
h_cell     = cell(1, cell_len);
record = zeros(1, cell_len);
for i = 1 : cell_len
    x         = org(1+(i-1)*L_frame+offset : 1+(i-1)*L_frame+offset+L_fagg-1);
    y         = rcv(1+(i-1)*L_frame : 1+(i-1)*L_frame+L_fagg-1);
    h         = get_CIR(x, y, a, b);
    h_cell{i} = h;
end

Gamma_cell3 = {};
N = length(h_cell{1});
tal = 1129;
for cnt = 1 : 10 : length(h_cell)-1
    disp(['h_cell no. ', int2str(cnt), ' total ', int2str(length(h_cell)-1)]);
    x = h_cell{cnt};
    xp = h_cell{cnt+1};
    Gamma = zeros(1, N);
    for k = 0 : N-1
        common = 0;
        for n = 0 : N-tal-1
            common = common + x(n+1)*exp(-1j*2*pi*k*n/N);
        end
        lower = common;
        for n = N-tal : N-1
            lower = lower + x(n+1)*exp(-1j*2*pi*k*n/N);
        end
        upper = common;
        for n = N-tal : N-1
            upper = upper + xp(n+1)*exp(-1j*2*pi*k*n/N);
        end
        Gamma(k+1) = upper / lower;
    end
    Gamma_cell3 = [Gamma_cell3, Gamma];
end

%%
fs         = 48000;       % sampling rante of a smartphone
T_frame    = 0.1;         % sampling period of frames  0.1 by default
T_fagg     = 0.2;         % compute CIR using T_fagg data    0.4 by default
L_frame    = fs*T_frame;  % length of a frame (CIR)
L_fagg     = fs*T_fagg;   % length of the aggregated data
music_code = 'orig0-';
echo_code  = 'XWT0_1';
clock_diff = 1;

% obtain coefficient k from .wav to .pcm
wav        = audioread('dataset/preamble.wav')'; %-1----clip of a .wav-----
pcm        = pcmread('dataset/preamble.pcm')'; %-1------clip of a .pcm-----
k          = max(pcm) / max(wav); %-1-------constant from .wav to .pcm-----

org        = audioread(['dataset/', music_code, '.wav'])';
rcv        = pcmread(['dataset/raw-data/', echo_code, '.pcm'])' / k; %-1---
% rcv        = audioread(['dataset/raw-data/', echo_code, '.wav'])'; %-2-----

% synchronization: sb_music -> music
[f1, l1]   = xcorr(rcv(1:150*fs), org(1:150*fs));        % correlate rcv with org
[~, idx1]  = max(f1);                                    % match point
rcv        = rcv(l1(idx1) + 1 : l1(idx1) + length(org)); % correlation

% 1.low pass filter
[b, a]     = butter(20, 15000/(fs/2));
% 2.find the index of the maxima in CIR
h_tmp      = zeros(1, L_fagg);
for i = 1 : L_frame : 1+50*L_frame
    x     = org(i+50 : i+L_fagg-1+50);
    y     = rcv(i : i+L_fagg-1);
    h_tmp = h_tmp + get_CIR(x, y, a, b);
end
maxima     = max(h_tmp);
dmax       = find(h_tmp == maxima);           % a CIR should start at 'dmax'
% dmax = 1;  %---------------------------------------------------------------
% 3.segmentation
offset     = 50;
cell_len   = ceil((length(org) - (L_fagg+offset) + 1) / (L_frame));
h_cell     = cell(1, cell_len);
record = zeros(1, cell_len);
for i = 1 : cell_len
    x         = org(1+(i-1)*L_frame+offset : 1+(i-1)*L_frame+offset+L_fagg-1);
    y         = rcv(1+(i-1)*L_frame : 1+(i-1)*L_frame+L_fagg-1);
    h         = get_CIR(x, y, a, b);
    h_cell{i} = h;
end

Gamma_cell2 = {};
N = length(h_cell{1});
tal = 1129;
for cnt = 1 : 10 : length(h_cell)-1
    disp(['h_cell no. ', int2str(cnt), ' total ', int2str(length(h_cell)-1)]);
    x = h_cell{cnt};
    xp = h_cell{cnt+1};
    Gamma = zeros(1, N);
    for k = 0 : N-1
        common = 0;
        for n = 0 : N-tal-1
            common = common + x(n+1)*exp(-1j*2*pi*k*n/N);
        end
        lower = common;
        for n = N-tal : N-1
            lower = lower + x(n+1)*exp(-1j*2*pi*k*n/N);
        end
        upper = common;
        for n = N-tal : N-1
            upper = upper + xp(n+1)*exp(-1j*2*pi*k*n/N);
        end
        Gamma(k+1) = upper / lower;
    end
    Gamma_cell2 = [Gamma_cell2, Gamma];
end