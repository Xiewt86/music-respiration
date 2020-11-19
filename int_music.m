function in = int_music (sig)
fs = 48000;
music = mapminmax(sig);
win = 0.1;
in = [];
for i = 1 : win*fs : length(music)-win*fs
    m = music(i : i+win*fs);
    M = abs(fft(m));
    in = [in, 20*log10(mean(M))];
end
end