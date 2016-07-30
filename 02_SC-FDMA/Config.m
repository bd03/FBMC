%% Config
%
% Burak Dayi
%
% Configuraton file in which the user may set number of subcarriers, 
% modulation type, number of symbols (number of frames), SNR value and
% number of trials.
%
% (For now it's only one-user case)
%
% In the configuration file user input should be controlled and if
% necessary, error messages should be raised.
%
% Last updated: 21-01-2015
c1 = clock;

%% General parameters
M = 256; % IFFT length
N = 12; % DFT length
mapping = 0; % 0: LFDMA 1:DFDMA
start = 15; % starting index for subcarrier mapping (1-M)
Q = ceil(M/N); % only for DFDMA mode bandwidth expansion factor

num_symbols = 50000 % number of symbols
modulation = 4 %4-, 16-, 64-, 128-, 256-QAM
num_samples = N %number of samples in a vector
bits_per_sample = log2(modulation) %num of bits carried by one sample
num_bits = num_symbols*num_samples*bits_per_sample % total number of bits transmitted
qam_sizes = [4 16 64 256];

cp_ratio = 0.25;%0.25;
cp_length = ceil(cp_ratio*M);

normalization = [2 10 42 170];

%---- Channel settings ----% 
% %rayleigh channel settings
% fading = 1; % set 0 for distortionless channel, set 1 for rayleigh channel
bw = 5e+6; % Transmission Bandwidth
% max_doppler_shift= 1; %max. doppler shift
% channel_profiles = ['EPA', 'EVA', 'ETU']; % Valid channel profile selections
% profile ='ETU'; %Channel profile
% use_matlab_channel = 0;
% if use_matlab_channel
%     [delay_a, pow_a] = LTE_channels2(profile,bw);
% %     ch_resp = rayleighchan(1/bw,max_doppler_shift,delay_a,pow_a); %channel model
%     ch_resp = stdchan(1/bw,max_doppler_shift,'itur3GVAx');
%     if max_doppler_shift>0
%         ch_resp.storeHistory = 1;
%     end
%     ch_resp.storePathGains =1;
%     ch_resp.ResetBeforeFiltering=0;
% else
%     ch_resp = LTE_channels(profile,bw);
% end
% 
% 
% %---- Channel estimation settings ----%
% % method
% valid_methods = ['IAM', 'IAM4']; % 'POP' may be added later.
% estimation_method = 'IAM';
% 
% % preamble will be defined during the simulation    
% 
% %---- Equalizer settings ----%
% eq_select = 2; % selection of equalizer type 1: one tap, 
% % 2: three tap w/ geometric interp, 3: three tap w/ linear interp
% % 4: no equalizer

% noise settings
noisy = 1; %set 1 for SNR values to affect channel, set 0 for noiseless channel
% SNR = 10; % SNR of the channel. noisy=1 to see the effects on channel
EbN0 = 0; % Eb/N0 ratio in dB.

noise_floor_dbHz = -165;
noise_floor_wattsHz = 10^(noise_floor_dbHz/10);
noise_pow = noise_floor_wattsHz*bw;

if (N>M)
    M
    N
    error('IFFT length M should be bigger than DFT length N.');
end

if ~(mapping ~= 0 | mapping ~= 1) 
    mapping
    error('Invalid mapping of symbols.');
end

if ~~mod(log2(M),1)
    M
    error('The number of subcarriers M should be a power of 2.');
end

if ~(ismember(modulation, qam_sizes))
    modulation
    %error('Only 4-QAM scheme is supported. Please set modulation = 4');
    error('Only 4-,16-,64-,128-,256-QAM schemes are supported. Define the number of constellation points (4,16,64,128,256) with modulation variable.')
end

%% Channel & simulation parameters
is_simulation = false;
SNR_array = 0:1:30; % Array of SNR value(s) in dB.
save('SNR_array.mat','SNR_array');
M_array = 2.^(10:10);
save('M_array.mat','M_array');
BER=zeros(length(M_array),length(qam_sizes),length(SNR_array));
num_trials = 100; % number of trials desired

if num_trials<1 || rem(num_trials,1)~=0
    num_trials
    error('Number of trials should be a positive integer.');
end

% if ~is_simulation
%     SNR = 1;
% end

ideal=0; %if channel is ideal

%% Additional configuration parameters
M_arr=M_array;
q_arr=qam_sizes;
s_arr=SNR_array;
if is_simulation
    conf=struct('M_val',M_arr,'mod_sch', q_arr,'SNR_val',s_arr);
    save(sprintf('CONF%d-%d-%d-%d-%d.mat', c1(1:5)),'conf');
end