#!/bin/env octave

n = 2^13
L = 2 * pi;
h = L/n;
x = (0:n-1)*h;
dt  = 0.0001;
mu = 0.001;
k = [0:n/2-1 0 -n/2+1:-1];
k1 = 1i*k;
k2 = k1.^2;
kf = 8 * pi;
				# uinit = cos(x.*x/3) / 2;
uinit = 1 - ((x - pi)/(pi)).^2;
uinit = center(uinit);
u0 = uinit;
E0 = sumsq(u0) * h;
q = 0;
figure('visible', 'off');
for t = 0:dt:2000
    E = sumsq(u0) * h;
    force = fft(u0);
    force(abs(k) > kf) = 0;
    force = 10 * force .* (E0 - E);
    u1hat = fft(u0) + 0.5*dt*(-k1.*fft(u0.^2)/2 + mu*k2.*fft(u0) + force);
    u1 = ifft(u1hat);
    u2hat = fft(u0) + dt*(-k1.*fft(u1.^2)/2 + mu*k2.*fft(u1) + force);
    u0 = ifft(u2hat);

    if mod(q, 1000) == 0
      plot(x, uinit, 'linewidth', 2, x, u0, 'linewidth', 5)
      axis([0,2*pi,-1, 1])
      title(sprintf('%12.1f', t))
      saveas(gcf(), sprintf('u.%09d.png', q))
      printf('%.16e\n', E);
      pause(0.01);
    end
    q += 1;
end
