%% Channel
%
% Burak Dayi
%
% disp('Channel')

y = tx_output;

% channel impulse response

if fading
    if use_matlab_channel
        y_filtered = filter(ch_resp,y);
    else
        y_filtered = filter(ch_resp,1,y);
    end    
else
    y_filtered = y;%filter(1,1,y);
end

% if ~noisy
%     y_ch = y_filtered;
% else
%     y_ch = awgn(y_filtered,SNR,'measured'); %SNR in dB ~ 10log(Ps/Pn)
% end

% mean_y_filtered=mean(abs(y_filtered.*y_filtered))

if ~noisy
    y_ch = y_filtered;
else
%     y_ch = awgn(y_filtered,SNR,'measured'); %SNR in dB ~ 10log(Ps/Pn)
	% sum_abs_yfilt=sum(abs(y_filtered.*y_filtered));
    amplitude = sqrt(num_bits*noise_floor_wattsHz*(10^(0.1*SNR))/sum(abs(y_filtered.*y_filtered)));
	% amplitude = sqrt(num_bits*noise_floor_wattsHz*bw*M/num_samples*(10^(0.1*SNR))/sum(abs(y_filtered.*y_filtered)));
    noise = randn(1,length(y_filtered))*sqrt(noise_pow/2)+j*randn(1,length(y_filtered))*sqrt(noise_pow/2);
    yy=amplitude*y_filtered;
    y_ch = yy + noise;
    mmm = find(qam_sizes==modulation);
    snrrr = find(EbN0_array==SNR);
    ebn0_check = 10*log10((sum(abs(yy.*yy))/num_bits)/noise_floor_wattsHz);
    % ebn0_check = 10*log10((sum(abs(yy.*yy))/num_bits)/mean(abs(noise.*noise)));
    ebno0_check_array(snrrr,mmm) = ebno0_check_array(snrrr,mmm)+(ebn0_check/num_trials);
end