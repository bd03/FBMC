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
% Created: 02-03-2014

% disp('Channel')

% channel impulse response
resp = [1]; % no impairments at the moment


if ideal
    y_ch = y;
else
    y_ch = awgn(y,SNR,'measured'); %SNR in dB ~ 10log(Ps/Pn)
end