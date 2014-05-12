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

if fading
    if use_matlab_channel
%         y_response = filter(ch_resp,ones(size(y)));
        y_filtered = filter(ch_resp,y);
%         y_filtered2 = y_response.*y;
    %     y_filtered = y_filtered2;
%         y_filtered3 = conv(y_response,y);
    else
        y_conv = conv(ch_resp,y);
        y_filtered = y_conv(length(ch_resp):length(ch_resp)+length(y)-1);  
    end    
else
    y_filtered = y;%filter(1,1,y);
end

if ~noisy
    y_ch = y_filtered;
else
    y_ch = awgn(y_filtered,SNR,'measured'); %SNR in dB ~ 10log(Ps/Pn)
end