clear all; close all; clc;

L = 30;
n = 512;

t2 = linspace(-L,L,n+1);
t = t2(1:n);
k = (2*pi)/(2*L)*[0:(n/2-1) -n/2:-1];

u = cos(2*t);

noise = 1;
ut = fft(u);
unt = ut + noise*(randn(1,n) + i*randn(1,n));
un = ifft(unt);

subplot(3,1,1), plot(t,u,'k'), hold on,
subplot(3,1,2),plot(t,abs(un),'k'),hold on,
subplot(3,1,3), plot(fftshift(k),abs(fftshift(unt))/max(abs(fftshift(unt))),'k'),
axis([-25 25 0 1]);