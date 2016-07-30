function [y_sum] = func_Channel(y,channel_params,oper_mode)
% function [y_ch] = func_Channel(y,fading,use_matlab_channel,profile,bw,noisy,SNR)
%% func_Channel
%
% Burak Dayi
%
% This function will perform channel impairments on the signal transmitted
% from the synthesis filter bank block.
%
% Dependencies: y - cell
% Output: y_ch - composite signal gone through channel response
%
% Created: 21-11-2014

% disp('Channel')

% channel impulse response

% first check if channel_params are all given.
try
    fading = channel_params.fading;
    use_matlab_channel = channel_params.use_matlab_channel;
    bw = channel_params.bw;
    profile = channel_params.profile;
    ch_resp = channel_params.ch_resp;
    noisy = channel_params.noisy;
    noise_floor_dbHz = channel_params.noise_floor_dbHz;
    noise_floor_wattsHz = channel_params.noise_floor_wattsHz;
    noise_pow = channel_params.noise_pow;
catch
    % missing element
    error('There are elements missing in channel_params structure.');
end

num_users = length(y);

y_filtered=cell(num_users,1);
y_sum = zeros(1,length(y{1})); % length assumed to be same for all: K*M+(num_oqam_subsymbols-1)*M/2

if strcmp(oper_mode,'MAIN')
    % In main mode, fading and noise variables are in form of vectors,
    % which means users might independently be exposed to fading.
    for i=1:num_users
        if fading(i)            
            if use_matlab_channel
                y_filtered = filter(ch_resp{i},y{i});
            else
                y_filtered = filter(ch_resp{i},1,y{i});
            end

            y_sum = y_sum + y_filtered;
        else
            y_sum = y_sum + y{i};
        end
    end
elseif strcmp(oper_mode,'SIMULATION')
    % apply channel if fading and sum
    if fading
        for i=1:num_users
            if use_matlab_channel
                y_filtered{i} = filter(ch_resp,y{i});
            else
                y_filtered{i} = filter(ch_resp,1,y{i});
            end

            y_sum = y_sum + y_filtered{i};
        end 
    else
        for i=1:num_users
            y_sum = y_sum + y{i};
        end
    end
else
    error('Undefined operation mode.')
end

% apply noise on signal.
noise = randn(1,length(y{1}))*sqrt(noise_pow/2)+j*randn(1,length(y{1}))*sqrt(noise_pow/2);
y_sum = y_sum + noise;

% noise = randn(1,K*M+(num_oqam_subsymbols-1)*M/2)*sqrt(noise_pow/2)+j*randn(1,K*M+(num_oqam_subsymbols-1)*M/2)*sqrt(noise_pow/2);



% if ~noisy
%     y_ch = y_sum;
% else
% %     y_ch = awgn(y_sum,SNR,'measured'); %SNR in dB ~ 10log(Ps/Pn)
%     % noise = randn(1,K*M+(num_oqam_subsymbols-1)*M/2)*sqrt(noise_pow/2)+j*randn(1,K*M+(num_oqam_subsymbols-1)*M/2)*sqrt(noise_pow/2);
%     y_ch = y_sum + noise;
%     mmm = find(qam_sizes==modulation);
%     snrrr = find(EbN0_array==SNR);
%     % ebn0_check = 10*log10((sum(abs(y{1}.*y{1}))/num_bits)/noise_floor_wattsHz);
%     % ebn0_check = 10*log10((sum(abs(yy.*yy))/num_bits)/mean(abs(noise.*noise)));
%     % ebno0_check_array(snrrr,mmm) = ebno0_check_array(snrrr,mmm)+(ebn0_check/num_trials);
% end