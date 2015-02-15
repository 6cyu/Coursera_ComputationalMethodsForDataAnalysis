%% Part 1
clear all; close all; clc;

load handel;
v = y'/2;
% plot((0:length(v)-1)/Fs,v);
% xlabel('Time [sec]');
% ylabel('Amplitude');
% title('Signal of Interest, v(n)');
p8 = audioplayer(v,Fs);
playblocking(p8);

L = (length(v)-1)/Fs;
n = length(v);
t2 = linspace(0,(length(v)-1)/Fs,n+1);
t = t2(1:n);
k = (1/L)*[0:n/2 -n/2:-1];
ks = fftshift(k);

vt = fft(v);

% width = [10 1 0.2 0.005];
width = [0.05 0.75 0.1 0.25];
tslide = 0:0.1:length(v)/Fs;
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
imagesc(tslide,ks,dum.'), shading interp,
ylim([-2000,2000]),
colormap(hot),colorbar;

dum(:,:) = spec(2,:,:);
imagesc(tslide,ks,dum.'), shading interp,
ylim([-2000,2000]),
colormap(hot),colorbar;

dum(:,:) = spec(3,:,:);
imagesc(tslide,ks,dum.'), shading interp,
ylim([-2000,2000]),
colormap(hot),colorbar;

dum(:,:) = spec(4,:,:);
imagesc(tslide,ks,dum.'), shading interp,
ylim([-2000,2000]),
colormap(hot),colorbar;

%% Part 2 - piano
clear all; close all; clc;
[y_piano2,Fs_piano] = audioread('music1.wav'); 
y_piano = y_piano2';
% Fs_piano = length(y_piano)/tr_piano;
n_piano = length(y_piano);
tr_piano = n_piano /Fs_piano;
t2_piano = linspace(0,n_piano,n_piano+1);
t_piano = t2_piano(1:n_piano)/Fs_piano;

k_piano = (1/tr_piano)*[0:n_piano/2-1 -n_piano/2:-1];
ks_piano = fftshift(k_piano);

% yt_piano = fft(y_piano);
plot(1:n_piano,y_piano);
xlabel('Time [sec]'); ylabel('Amplitude');
title('Mary had a little lamb (piano)'); drawnow;
% p8_piano = audioplayer(y_piano,Fs_piano); playblocking(p8_piano);

w_piano = [0.25 0.5 0.75 1];
trans = 0.1;
slide_piano = 0:trans:tr_piano;
spec_piano = zeros(length(w_piano),length(slide_piano),n_piano);

for i=1:length(w_piano)
    for j=1:length(slide_piano)
%         f_piano = exp(-w_piano(i)*(t_piano-slide_piano(j)).^2);
        f_piano = heaviside(t_piano - slide_piano(j)) - ...
            heaviside(t_piano - w_piano(i) - slide_piano(j));
        
        yf_piano = y_piano.*f_piano;
        
        yft_piano = fft(yf_piano);
        
        spec_piano(i,j,:) = abs(fftshift(yft_piano));
        
%         figure(3)
%         subplot(3,1,1), plot(t_piano,y_piano,t_piano,f_piano,'r');
%         figure(3)
%         subplot(3,1,2), plot(t_piano, yf_piano);
%         figure(3)
%         subplot(3,1,3), plot(ks_piano,...
%             abs(fftshift(yft_piano))/max(abs(fftshift(yft_piano))));
        
    end
end

%% Plotting piano spectrogram
% figure(4)
%%
dum(:,:) = spec_piano(1,:,:);
figure(4), 
% subplot(2,2,1), 
% pcolor(slide_piano,ks_piano,dum.'), shading interp,
imagesc(slide_piano,ks_piano,dum.'), shading interp,
ylim([-500,500]),
colormap(hot),colorbar;
% saveas(gcf,'H_S01W025','tif')
% exportfig('H_S01W025.png',...
%     'Width',10,...
%     'color','rgb'...
%     );


%%
dum(:,:) = spec_piano(2,:,:);
figure(5), 
% subplot(2,2,2),
% pcolor(slide_piano,ks_piano,dum.'), shading interp,
imagesc(slide_piano,ks_piano,dum.'), shading interp,

colormap(hot),colorbar;

%%
dum(:,:) = spec_piano(3,:,:);
figure(6), 
% subplot(2,2,3),
% pcolor(slide_piano,ks_piano,dum.'), shading interp,
imagesc(slide_piano,ks_piano,dum.'), shading interp,
colormap(hot),colorbar;
%%
dum(:,:) = spec_piano(4,:,:);
figure(7), 
% subplot(2,2,4),
% pcolor(slide_piano,ks_piano,dum.'), shading interp,
imagesc(slide_piano,ks_piano,dum.'), shading interp,
ylim([-500,500]),
colormap(hot),colorbar;
%         
% 
% exportfig('H_S01W025.png',...
%     'width',3.7,...
%     'color','rgb'...
%     );


%% Part 2 - recorder
clear all; close all; clc;
[y_rec2,Fs_rec] = audioread('music2.wav'); 
y_rec = y_rec2';
% Fs_rec = length(y_rec)/tr_rec;
n_rec = length(y_rec);
tr_rec = n_rec /Fs_rec;
t2_rec = linspace(0,n_rec,n_rec+1);
t_rec = t2_rec(1:n_rec)/Fs_rec;

k_rec = (1/tr_rec)*[0:n_rec/2-1 -n_rec/2:-1];
ks_rec = fftshift(k_rec);

% yt_rec = fft(y_rec);
plot(1:n_rec,y_rec);
xlabel('Time [sec]'); ylabel('Amplitude');
title('Mary had a little lamb (rec)'); drawnow;
% p8_rec = audioplayer(y_rec,Fs_rec); playblocking(p8_rec);

w_rec = [0.25 0.5 0.75 1];
trans = 0.1;
slide_rec = 0:trans:tr_rec;
spec_rec = zeros(length(w_rec),length(slide_rec),n_rec);

for i=1:length(w_rec)
    for j=1:length(slide_rec)
        f_rec = heaviside(t_rec - slide_rec(j)) - ...
            heaviside(t_rec - w_rec(i) - slide_rec(j));
        
        yf_rec = y_rec.*f_rec;
        
        yft_rec = fft(yf_rec);
        
        spec_rec(i,j,:) = abs(fftshift(yft_rec));
        
%         figure(3)
%         subplot(3,1,1), plot(t_rec,y_rec,t_rec,f_rec,'r');
%         figure(3)
%         subplot(3,1,2), plot(t_rec, yf_rec);
%         figure(3)
%         subplot(3,1,3), plot(ks_rec,...
%             abs(fftshift(yft_rec))/max(abs(fftshift(yft_rec))));
        
    end
end

%% Plotting piano spectrogram
% figure(4)
%%
dum(:,:) = spec_rec(1,:,:);
figure(4), 
% subplot(2,2,1), 
% pcolor(slide_piano,ks_piano,dum.'), shading interp,
imagesc(slide_rec,ks_rec,dum.'), shading interp,
ylim([-1500,1500]),
colormap(hot),colorbar;
% saveas(gcf,'H_S01W025','tif')
% exportfig('H_S01W025.png',...
%     'Width',10,...
%     'color','rgb'...
%     );


%%
dum(:,:) = spec_rec(2,:,:);
figure(5), 
% subplot(2,2,2),
% pcolor(slide_piano,ks_piano,dum.'), shading interp,
imagesc(slide_rec,ks_rec,dum.'), shading interp,
ylim([-500,500]),
colormap(hot),colorbar;

%%
dum(:,:) = spec_rec(3,:,:);
figure(6), 
% subplot(2,2,3),
% pcolor(slide_piano,ks_piano,dum.'), shading interp,
imagesc(slide_rec,ks_rec,dum.'), shading interp,
ylim([-500,500]),
colormap(hot),colorbar;
%%
dum(:,:) = spec_rec(4,:,:);
figure(7), 
% subplot(2,2,4),
% pcolor(slide_piano,ks_piano,dum.'), shading interp,
imagesc(slide_rec,ks_rec,dum.'), shading interp,
ylim([-500,500]),
colormap(hot),colorbar;
%         
% 
% exportfig('H_S01W025.png',...
%     'width',3.7,...
%     'color','rgb'...
%     );
