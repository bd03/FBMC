%% Channel
%
% Burak Dayi
%
% This script will perform channel impairments on the signal transmitted
% from the synthesis filter bank block.
%
% Dependencies: y - composite signal output
% Output: y_ch - composite signal gone through channel response
%
% Last updated: 18-03-2014

% disp('Channel')

% channel impulse response
resp = [1]; % no impairments at the moment

%y_ch = y; ideal=1;
y_ch = awgn(y,SNR,'measured'); ideal=0; %SNR in dB ~ 10log(Ps/Pn)