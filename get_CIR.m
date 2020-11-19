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