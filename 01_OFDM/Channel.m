%% Channel
%
% Burak Dayi
%
% disp('Channel')

y = tx_output;

% channel impulse response
resp = [1]; % no impairments at the moment

%y_ch = y; ideal=1;
y_ch = awgn(y,SNR,'measured'); ideal=0; %SNR in dB ~ 10log(Ps/Pn)