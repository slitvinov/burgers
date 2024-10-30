#!/bin/env octave

randn('state', 1);
n = 2^13;
L = 2 * pi;
h = L/n;
x = (0:n-1)*h;
dt  = 0.0001;
mu = 0.001;
k = [0:n/2 -n/2+1:-1];
k1 = 1i*k; k1(n/2+1) = 0;
k2 = -k.^2;

U = k .^2 .* exp(-k.^2 / L);
a =  fft(randn(1, n));
uinit = ifft(U .* a);
uinit = 0.8 * uinit / max(abs(uinit));
u0 = uinit;
E0 = sumsq(u0) * h;
En = zeros(size(k));
q = 0;
figure('visible', 'off');
for t = 0:dt:200000
    E = sumsq(u0) * h;
    force = fft(u0);
    force(~(abs(k) <= 4)) = 0;
    force .*= (E0 - E) / dt;
    rhs = @(u) (-k1.*fft(u0.^2)/2 + mu*k2.*fft(u0) + force);

    u1hat = fft(u0) + 0.5*dt*rhs(u0);
    u1hat(abs(k) > n/2 * 2/3) = 0;

    u1 = ifft(u1hat);
    u2hat = fft(u0) + dt*rhs(u1);
    u2hat(abs(k) > n/2 * 2/3) = 0;

    u0 = ifft(u2hat);
    q += 1;
    En = En * (q - 1) / q  + abs(u2hat).^2/2*h / q;

    if mod(q, 1000) == 0
      uf = u2hat;
      uf(abs(k) >= 70) = 0;
      uf = ifft(uf);
      force = ifft(force);
      force /= max(abs(force));
      plot(x, u0, 'color', 'black', 'linewidth', 8,
	   x, uf, 'color', 'black', 'linewidth', 2,
	   x, force, 'color', 'black', 'linewidth', 2)
      axis([0 2*pi -1 1])
      title(sprintf('%12.2f', t))
      saveas(gcf(), sprintf('u.%09d.png', q))

      % loglog(k(k > 0), h * En(k > 0) / 2, 'b.-');
      % title(sprintf('%12.2f', t))
      % saveas(gcf(), sprintf('e.%09d.png', q))

      printf('%.16e %.16e\n', E, (E - E0) / (mu^3));
    end
end
