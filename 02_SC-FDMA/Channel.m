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
%
% Created on 30-01-2015

% if fading
%     if use_matlab_channel
%         y_filtered = filter(ch_resp,y);
%     else
%         y_filtered = filter(ch_resp,1,y);
%     end    
% else
%     y_filtered = y;  % filter(1,1,y);
% end
% 
if ~noisy
    % y_ch = y_filtered;
    y_ch = tx_output;
else
    % y_ch = awgn(y_filtered,SNR,'measured'); %SNR in dB ~ 10log(Ps/Pn)
    % y_ch = awgn(tx_output,SNR,'measured'); %SNR in dB ~ 10log(Ps/Pn)
    noise = randn(1,(1+cp_ratio)*M*num_symbols)*sqrt(noise_pow/2)+j*randn(1,(1+cp_ratio)*M*num_symbols)*sqrt(noise_pow/2);
    y_ch = noise + tx_output;
end

% y_ch = tx_output;