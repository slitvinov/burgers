#!/bin/env octave

graphics_toolkit('gnuplot');
randn('state', 1);
n = 2^13;
L = 2 * pi;
h = L/n;
x = (0:n-1)*h;
dt  = 0.0001;
mu = 0.001;
k = [0:n/2-1 -n/2:-1];
k1 = 1i*k;
k2 = k1.^2;

U = k .^2 .* exp(-k.^2 / L);
a =  fft(randn(1, n));
uinit = ifft(U .* a);
uinit = 0.8 * uinit / max(abs(uinit));
u0 = uinit;
E0 = sumsq(u0) * h;
q = 0;
figure('visible', 'off');
for t = 0:dt:2000
    E = sumsq(u0) * h;
    force = fft(u0);
    force(abs(k) >= 3) = 0;
    force .*= (E0 - E) / dt;
    rhs = @(u) (-k1.*fft(u0.^2)/2 + mu*k2.*fft(u0) + force);
    u1hat = fft(u0) + 0.5*dt*rhs(u0);
    u1 = ifft(u1hat);
    u2hat = fft(u0) + dt*rhs(u1);
    u0 = ifft(u2hat);
    if mod(q, 5000) == 0
      uf = u2hat;
      uf(abs(k) >= 70) = 0;
      uf = ifft(uf);
      force = ifft(force);
      force /= max(abs(force));
      plot(x, u0, 'color', 'black', 'linewidth', 8,
	   x, uf, 'color', 'black', 'linewidth', 2,
           x, force, 'color', 'black', 'linewidth', 2)
      axis([0,2*pi,-1, 1])
      title(sprintf('%12.2f', t))
      saveas(gcf(), sprintf('u.%09d.png', q))
      printf('%.16e\n', E);
    end
    q += 1;
end
