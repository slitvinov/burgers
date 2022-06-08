n = 2^11;
L = 2 * pi;
h = L/n;
x = (0:n-1)*h;
k = [0:n/2-1 0 -n/2+1:-1];
U = 1/sqrt(3) * (2 * pi)^(1/4) * k .^2 .* exp(-k.^2 / 4);
a =  fft(randn(1, n));
uinit = ifft(U .* a) * 8;
plot(x, uinit)
axis([0,2*pi,-1, 1])
pause
