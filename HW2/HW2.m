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
[y_piano,Fs_piano] = audioread('music1.wav'); 
% Fs_piano = length(y_piano)/tr_piano;
tr_piano = length(y_piano)/Fs_piano;
n_piano = 
plot((1:length(y_piano))/Fs_piano,y_piano);
xlabel('Time [sec]'); ylabel('Amplitude');
title('Mary had a little lamb (piano)'); drawnow;
p8_piano = audioplayer(y_piano,Fs_piano); playblocking(p8_piano);




%% Part 2 - recorder

figure(2)
[y_rec,Fs_rec] = audioread('music2.wav'); 
% Fs_rec = length(y_rec)/tr_rec;
tr_rec = length(y_rec)/Fs_rec;
plot((1:length(y_rec))/Fs_rec,y_rec);
xlabel('Time [sec]'); ylabel('Amplitude');
title('Mary had a little lamb (recorder)');  
p8_rec = audioplayer(y_rec,Fs_rec); playblocking(p8_rec);
