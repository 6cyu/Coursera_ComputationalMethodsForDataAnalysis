clear all; close all; clc;

L = 10;
n = 2048;
t2 = linspace(0,L,n+1); t = t2(1:n);
k = (2*pi/L)*[0:n/2-1 -n/2:-1];
ks = fftshift(k);

E = (3*sin(2*t) + 0.5*tanh(0.5*(t-3))+...
    0.28*exp(-(t-4).^2)...
    + 1.5*sin(5*t)+4*cos(3*(t-6).^2))/10 + ...
    (t/20).^3;

Et = fft(E);

width = [10 1 0.2];
for j=1:3
    f = exp(-width(j)*(t-4).^2);
    subplot(3,1,j), plot(t,E,'k',t,f,'m');
    
end

width = 0.1;
slide = 0:0.1:10;
spec = [];
for j=1:length(slide)
    f = exp(-width*(t-slide(j)).^2);
    
    
    % signal filtered
    Ef = E.*f;
    Eft = fft(Ef);
    
    subplot(3,1,1), plot(t,E,'k',t,f,'m');
    subplot(3,1,2), plot(t,Ef,'k');
    subplot(3,1,3), plot(ks,abs(fftshift(Eft))/...
                        max(abs(fftshift(Eft))));
    axis([-60 60 0 1]);
    drawnow;
    pause(0.1);
    
    spec = [spec; abs(fftshift(Eft))];
    
end

figure,pcolor(slide,ks,spec.'), shading interp;
set(gca,'Ylim',[-60 60], 'Fontsize', [14]);
colormap(hot);
xlabel('t');
ylabel('omega');
colorbar;


subplot(3,1,1), plot(t,E1,t,f1);
subplot(3,1,1), plot(t,E2,t,f2);
subplot(3,1,1), plot(t,E3,t,f3);
subplot(2,1,2), plot(ks,...
    abs(fftshift(Et))/(max(abs(fftshift(Et)))));



