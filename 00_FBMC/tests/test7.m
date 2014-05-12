%test7.m

% test to get freq coefficients of an LTE channel

clear all
close all
clc
bw = 5e+6;
nsym=1e+3;
max_doppler_shift= 1e-5;
channel_profiles = ['EPA', 'EVA', 'ETU']; % Valid channel profile selections
profile= 'ETU';
ch_resp = LTE_channels (profile,bw);
x=zeros(1,nsym);
x(1)=1;
y= conv(ch_resp,x);
% y_filtered = y(length(ch_resp):length(ch_resp)+length(y)-1);  
clear x;
fft_y=fft(y);
z=zeros(1,length(y));
z(1)=1;
fft_t=fft(conv(z,y));
clear y;

% fft_y=[fft_y(floor(length(fft_y)/2)+1:end) fft_y(1:floor(length(fft_y)/2))];

%plotting
figure;
plot(abs(fft(ch_resp,length(fft_y))))
% axis([0 nsym -40 10])
% axis([0 nsym 0.96 1.04])
title('ch_resp')
grid on
figure;
plot(abs(fft_y))
% axis([0 nsym -40 10])
% axis([0 length(fft_y) 0.96 1.04])
grid on
title('fft_y')
figure;
plot(abs(fft_t))
% axis([0 nsym -40 10])
% axis([0 length(fft_y) 0.96 1.04])
grid on
title('fft_t')

mser=sum(abs((fft_y-fft(ch_resp,length(fft_y))).^2))/size(fft_y,1)