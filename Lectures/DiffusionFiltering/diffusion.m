clear all; close all; clc;

A = imread('lena','png');
Abw = rgb2gray(A);

A2 = double(Abw);
[nx,ny] = size(A2);

An = A2 + 50*randn(nx,ny);

x = linspace(0,1,nx); % normalize frequency domain
y = linspace(0,1,ny);
dx = x(2)-x(1);
dy = y(2)-y(1);

onex = ones(n,1); oney = ones(ny,1);
Dx = spdiags([onex -2*onex onex],[-1 0 1], nx,nx)/dx^2;
Dy = spdiags([oney -2*oney oney],[-1 0 1], ny,ny)/dy^2;

Ix = eye(nx); Iy = eye(ny);
L = kron(Iy,Dx) + kron(Dy,Ix);

spy(L);

tspan = [ 0 0.002 0004 0006];
An2 = reshape(An,nx*ny,1);

D = 1;
[t, uso] = ode45('zoo_rhs',tspan,An2,[],L,D);

for j)1>length*t(
    Atemp ) uint8*reshape*uso*j(