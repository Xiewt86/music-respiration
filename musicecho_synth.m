%------------------------synthesis signal------------------------------%

fs = 44100;
load('test1frame.mat');
echo = pcmread('th3.pcm');
receive_OFDM_cont;
h = hh{1}(1 : length(hh{1})/45);

music_wav = audioread('music.wav')';
music_pcm = pcmread('music.pcm')';
music_echo = pcmread('sb8.pcm');

% coefficient: from .wav to .pcm
k = max(music_pcm) / max(music_wav);

% divided by k
h = h / k;
music_echo = music_echo / k;

% synthetic signal
synth = conv(music_wav, h);

% synchronization: music_echo -> music
[f0, lag0] = xcorr(music_echo, music_wav);
[~, idx0] = max(f0);
music_echo = music_echo( lag0(idx0)+1 : lag0(idx0)+length(music_wav) );

% synchronization: synth -> music_echo
[f1, lag1] = xcorr(synth, music_echo);
[~, idx1] = max(f1);
synth = synth( lag1(idx1)+1 : lag1(idx1)+length(music_echo) );

% LPF
[b, a] = butter(8, 200/(fs/2), 'high');
% music_echo = filter(b, a, music_echo);
% synth = filter(b, a, synth);

synth0 = energy_adjust(music_echo, synth);
figure; plot(synth0);
hold on; plot(music_echo);
legend('synth', 'echo')

figure; plot(music_echo);
hold on; plot(synth0);
legend('echo', 'synth')

% no_music = synth(1:end)';
% sb_music = music_echo(1:end)';

