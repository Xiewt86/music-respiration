load('Gammas.mat');
snr4 = 0;
for i = 1 : length(Gamma_cell4)
    G = Gamma_cell4{i};
    g = abs(ifft(G(2:end)));
    snr4 = snr4 + 1 / (1-g(1)^2/sum(g.^2));
end
snr4 = snr4 / length(Gamma_cell4);

snr3 = 0;
for i = 1 : length(Gamma_cell3)
    G = Gamma_cell3{i};
    g = abs(ifft(G(2:end)));
    snr3 = snr3 + 1 / (1-g(1)^2/sum(g.^2));
end
snr3 = snr3 / length(Gamma_cell3);
    
snr2 = 0;
for i = 1 : length(Gamma_cell2)
    G = Gamma_cell2{i};
    g = abs(ifft(G(2:end)));
    snr2 = snr2 + 1 / (1-g(1)^2/sum(g.^2));
end
snr2 = snr2 / length(Gamma_cell2);
