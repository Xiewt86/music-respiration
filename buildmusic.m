fs = 48000;
p = audioread('preamble(2).wav');

musict = p';
musict = [musict, zeros(1, 0.2*fs)];
music = [];

for i = 1 : 6
    name = ['a', mat2str(i), '(2).wav'];
    m = audioread(name);
    musict = [musict, m'];
    music = [music, m'];
end

audiowrite('orig1.wav', music, fs)
audiowrite('testing1.wav', musict, fs)

musict = p';
musict = [musict, zeros(1, 0.2*fs)];
music = [];

for i = 1 : 6
    name = ['b', mat2str(i), '(2).wav'];
    m = audioread(name);
    musict = [musict, m'];
    music = [music, m'];
end

audiowrite('orig2.wav', music, fs)
audiowrite('testing2.wav', musict, fs)

