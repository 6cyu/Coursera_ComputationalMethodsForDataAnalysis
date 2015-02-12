clear all; close all; clc;

load Testdata;

L = 15; %spatial domain
n = 64; %Fourier domain
x2 = linspace(-L,L,n+1);
x = x2(1:n);
y = x;
z = x;
k = (2*pi/L) * [0:(n/2-1) -n/2:-1];
ks = fftshift(k);

[X,Y,Z] = meshgrid(x,y,z);
[Kx,Ky,Kz] = meshgrid(ks,ks,ks);

% UntAvg = complex(zeros(n,n,n));
Un = zeros(n,n,n);
Unt = zeros(n,n,n);
UntAvg = zeros(n,n,n);

for j=1:20
    Un = reshape(Undata(j,:),n,n,n);        
    Unt = fftn(Un);
    UntAvg = UntAvg + Unt;
end

UntAvg = UntAvg./20;
[maxAvg, maxIndexAvg] = max(abs((UntAvg(:))));
[yIndexAvg, xIndexAvg, zIndexAvg] = ind2sub(size(UntAvg), maxIndexAvg);

isovalue = 200;

KxC = k(xIndexAvg);
KyC = k(yIndexAvg);
KzC = k(zIndexAvg);

isosurface(Kx,Ky,Kz,abs(fftshift(UntAvg)),200);

% Gaussian filter
sigma = 0.35;
% GF = exp(-((Kx-KxC).^2 + (Ky-KyC).^2 + (Kz-KzC).^2)/sigma);
fx = exp(-((k-KxC).^2)/sigma);
fy = exp(-((k-KyC).^2)/sigma);
fz = exp(-((k-KzC).^2)/sigma);
[FX, FY, FZ] = meshgrid(fx, fy, fz);

Untf20 = ((Unt.*FX).*FY).*FZ;
% Untf20 = Unt.*GF;
UInv20 = ifftn((Untf20));

for j=0:10
    figure,isosurface(X,Y,Z,abs((UInv20)),j*0.01);
end
isosurface(X,Y,Z,abs((UInv20)),0.01);
axis([-20 20 -20 20 -20 20]),grid on, drawnow;

[maxA, maxIndex] = max(abs((UInv20(:))));
[yIndex, xIndex, zIndex] = ind2sub(size(UInv20), maxIndex);
x20 = x(xIndex);
y20 = x(yIndex);
z20 = x(zIndex);

close all;

for j=1:20
    Un = reshape(Undata(j,:),n,n,n);        
    Unt = fftn(Un);
    Untf = ((Unt.*FX).*FY).*FZ;
%     Untf = Unt.*GF;
    Unf = ifftn(Untf);
    figure,isosurface(X,Y,Z,abs((Unf)),0.4);
    axis([-20 20 -20 20 -20 20]),grid on, drawnow;
end