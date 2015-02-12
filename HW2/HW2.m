%% Part 1
clear all; close all; clc;

load handel;
v = y'/2;
% plot((0:length(v)-1)/Fs,v);
% xlabel('Time [sec]');
% ylabel('Amplitude');
% title('Signal of Interest, v(n)');
% p8 = audioplayer(v,Fs);
% playblocking(p8);

L = (length(v)-1)/Fs;
n = length(v);
t2 = linspace(0,(length(v)-1)/Fs,n+1);
t = t2(1:n);
k = (2*pi/L)*[0:n/2 -n/2:-1];
ks = fftshift(k);

vt = fft(v);

% width = [10 1 0.2 0.005];
width = [0.25 0.5 0.75 1];
tslide = 0:0.01:length(v)/Fs;
spec = zeros(length(width),length(tslide),n);
for i=1:length(width)
    for j=1:length(tslide)
%         f = exp(-width(i)*(t-tslide(j)).^2);
        f = heaviside(t-tslide(j)) - heaviside(t-width(i)-tslide(j));
        vf = v.*f;
        vft = fft(vf);
        spec(i,j,:) = abs(fftshift(vft));
        
%         figure(1),
%         subplot(3,1,1),plot(t,v,t,f,'r');
%         figure(1),
%         subplot(3,1,2),plot(t,vf);
%         figure(1),
%         subplot(3,1,3),plot(ks,abs(fftshift(vft))/max(abs(fftshift(vft))));
%         pause(0.1);  
        
    end
end

%% spectrogram plot
dum(:,:) = spec(1,:,:);
figure(2), subplot(2,2,1), 
pcolor(tslide,ks,dum.'), shading interp,
colormap(hot),colorbar;
dum(:,:) = spec(2,:,:);
figure(2), subplot(2,2,2),
pcolor(tslide,ks,dum.'), shading interp,
colormap(hot),colorbar;
dum(:,:) = spec(3,:,:);
figure(2), subplot(2,2,3),
pcolor(tslide,ks,dum.'), shading interp,
colormap(hot),colorbar;
dum(:,:) = spec(4,:,:);
figure(2), subplot(2,2,4),
pcolor(tslide,ks,dum.'), shading interp,
colormap(hot),colorbar;

%% Part 2 - piano
[y_piano2,Fs_piano] = audioread('music1.wav'); 
y_piano = y_piano2';
% Fs_piano = length(y_piano)/tr_piano;
n_piano = length(y_piano);
tr_piano = n_piano /Fs_piano;
t2_piano = linspace(0,n_piano,n_piano+1);
t_piano = t2_piano(1:n_piano);

k_piano = (2*pi/tr_piano)*[0:n_piano/2-1 -n_piano/2:-1];
ks_piano = fftshift(k_piano);

yt_piano = fft(y_piano);

plot(1:n_piano,y_piano);
xlabel('Time [sec]'); ylabel('Amplitude');
title('Mary had a little lamb (piano)'); drawnow;
% p8_piano = audioplayer(y_piano,Fs_piano); playblocking(p8_piano);

w_piano = [0.25 0.5 0.75 1];
slide_piano = 0:0.1:tr_piano;
spec_piano = zeros(length(w_piano),length(slide_piano),n_piano);

for i=1:length(w_piano)
    for j=1:length(slide_piano)
        f_piano = heaviside(t_piano - slide_piano(j)) - ...
            heaviside(t_piano - w_piano(i) - slide_piano(j));
        
        yf_piano = y_piano.*f_piano;
        
        yft_piano = fft(yf_piano);
        
        spec_piano(i,j,:) = yft_piano;
        
        figure(3)
        subplot(3,1,1), plot(t_piano,y_piano,t_piano,f_piano,'r');
        figure(3)
        subplot(3,1,2), plot(t_piano, yf_piano);
        figure(3)
        subplot(3,1,3), plot(ks_piano,...
            abs(fftshift(yft_piano))/abs(fftshift(yft_piano)));
        
    end
end

%% Plotting piano spectrogram
figure(4)
dum(:,:) = spec_piano(1,:,:);
figure(4), subplot(2,2,1), 
pcolor(slide_piano,ks_piano,dum.'), shading interp,
colormap(hot),colorbar;
dum(:,:) = spec_piano(2,:,:);
figure(4), subplot(2,2,2),
pcolor(slide_piano,ks_piano,dum.'), shading interp,
colormap(hot),colorbar;
dum(:,:) = spec_piano(3,:,:);
figure(4), subplot(2,2,3),
pcolor(slide_piano,ks_piano,dum.'), shading interp,
colormap(hot),colorbar;
dum(:,:) = spec_piano(4,:,:);
figure(4), subplot(2,2,4),
pcolor(slide_piano,ks_piano,dum.'), shading interp,
colormap(hot),colorbar;
        




%% Part 2 - recorder

figure(2)
[y_rec,Fs_rec] = audioread('music2.wav'); 
% Fs_rec = length(y_rec)/tr_rec;
n_rec = length(y_rec);
tr_rec = n_rec/Fs_rec;
plot((1:n_rec)/Fs_rec,y_rec);
xlabel('Time [sec]'); ylabel('Amplitude');
title('Mary had a little lamb (recorder)');  
p8_rec = audioplayer(y_rec,Fs_rec); playblocking(p8_rec);
