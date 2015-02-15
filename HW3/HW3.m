%% Derek1 RGB
clear all; close all; clc;

D= imread('derek1','jpg');
Dbw = rgb2gray(D);

D2 = im2double(Dbw);
[nx,ny] = size(D2);

D2t = fft2(D2);
D2ts = fftshift(D2t);

%% constructing gaussian filter
kx = 1:ny;
ky = 1:nx;
[Kx,Ky] = meshgrid(kx,ky);
wG = 1e-5;
FG = exp(-wG*((Kx-ceil(ny/2)).^2 - (Ky-ceil(nx/2)).^2));
% pcolor(f),colormap(hot),shading interp;

%% constructing Shannon filter
wS = 75;
FS = zeros(size(D2));
FS(ceil(nx/2)-wS:1:ceil(nx/2)+wS, ceil(ny/2)-...
    wS:1:ceil(ny/2)+wS)...
    = ones(2*wS+1,2*wS+1);

% pcolor(FS),colormap(hot),shading interp;


%% filtering with gaussian filter
D2tsf = D2ts.*FG;
D2tf = ifftshift(D2tsf);
D2f = (ifft2((D2tf)));


%% filtering with shannon filter
D2tsf = D2ts.*FS;
D2tf = ifftshift(D2tsf);
D2f = ifft2(D2tf);


%% displaying image
subplot(2,2,1), imshow(D2,[]), colormap(gray);
subplot(2,2,2), pcolor(log(abs(D2ts))), shading interp,...
    colormap(hot), set(gca,'Xtick',[],'Ytick',[]);
subplot(2,2,3), pcolor(log(abs(D2tsf))), shading interp,...
    colormap(hot), set(gca,'Xtick',[],'Ytick',[]);
subplot(2,2,4), imshow((D2f),[]), colormap(gray);

%% converting to rgb
rgbD2f = repmat(D2f,[1 1 3]);
rgbD2f = cat(3,rgbD2f,rgbD2f,rgbD2f);
imshow(rgbD);
gray2rgb

%%



x = linspace(0,1,nx); % normalize frequency domain
y = linspace(0,1,ny);
dx = x(2)-x(1);
dy = y(2)-y(1);

onex = ones(ny,1); oney = ones(ny,1);
Dx = spdiags([onex -2*onex onex],[-1 0 1], nx,nx)/dx^2;
Dy = spdiags([oney -2*oney oney],[-1 0 1], ny,ny)/dy^2;

Ix = eye(nx); Iy = eye(ny);
L = kron(Iy,Dx) + kron(Dy,Ix);

spy(L);

tspan = [ 0 0.002 0004 0006];
An2 = reshape(An,nx*ny,1);

D = 1;
[t, uso] = ode45('zoo_rhs',tspan,An2,[],L,D);
% 
% for j)1>length*t(
%     Atemp ) uint8*reshape*uso*j(