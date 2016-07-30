function [h] = func_Prototype_filter(M,K,lp)
%% Prototype_filter
%
% Burak Dayi
%
% This script will compute the prototype filter coefficients and then plot 
% the prototype filter in frequency and time domain.
%
% Dependencies: lp, K, M, delay, P
% Output: h - prototype filter coefficients in time domain
%
% Created: 21-11-2014

% disp('Prototype filter design')

% computed again because we have to modify lp for each M value in
% FBMC_Simulation
%
% ---- Prototype filter frequency coefficients----%
%---------------------------------
% This part should not be altered!
%---------------------------------
% K will select the row, the last column is background noise power in dB 
% (that might be come in handy in future) for reference
P=[ zeros(1,5);1 sqrt(2)/2 0 0 -35; 1 .911438 .411438 0 -44; 1 .97195983 sqrt(2)/2 .23514695 -65];

delay = K*M+1-lp; %delay requirement 

h = zeros(1,lp); %coefficient array created with filter length
for m=1:lp
    h(m) = P(K,1);
    for k = 1:(K-1)
        h(m) = h(m) + 2*(((-1)^k)*P(K,k+1)*cos((2*pi*k*m)/(K*M)));
    end
end

%normalization of coefficients
denom = P(K,1)+2*sum(P(K,2:4));
h=h/denom; %normalization

%delay implementation as on p.19 of deliverable 5.1
if delay>1
    h = [0 h 0];
    delay = delay - 1;
end

% figure(11);
% subplot(311);
% plot(h);
% title(sprintf('Time response of PHYDYAS Filter with K=%d M=%d lp=%d', K,M,lp));
% grid on;
% subplot(312);
% plot(abs(fft(h)));
% title(sprintf('Freq response of PHYDYAS Filter with K=%d M=%d lp=%d', K,M,lp));
% ylabel('Magnitude');
% grid on;
% subplot(313);
% semilogy(abs(fft(h)));
% title(sprintf('Freq magnitude response of PHYDYAS Filter with K=%d M=%d lp=%d', K,M,lp));
% % ylabel('Log Magnitude');
% grid on;
% close all;