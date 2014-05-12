%test6.m

% test to get freq coefficients of a built-in matlab channel

clear all
close all
clc
bw = 5e+6;
nsym=1e+3;
max_doppler_shift= 1e-5;
ch_resp = stdchan(1/bw,max_doppler_shift,'itur3GVAx');
ch_resp.StoreHistory=1;
x=zeros(1,nsym);
x(1)=1;
y=filter(ch_resp,x);
% y=[y(floor(length(y)/2+1:end)) y(1:floor(length(y)/2))] ;
clear x;
fft_y=fft(y,length(y)); %+size(ch_resp.PathGains,2)
clear y;
fft_y=[fft_y(floor(length(fft_y)/2)+1:end) fft_y(1:floor(length(fft_y)/2))];

%plotting
plot(ch_resp)
plot(20*log10(abs(fft_y)))
axis([0 nsym -40 10])
grid on