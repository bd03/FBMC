%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Config_Multiuser.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created: 05-02-2014

is_simulation = false;

%--------General Parameters--------%
M = 256; % IFFT length
N = 12; % DFT length for each user
num_symbols = 500; % number of symbols
modulation = 64; %4-, 16-, 64-, 128-, 256-QAM
allocated_subchannels = N; %number of samples in a vector
bits_per_sample = log2(modulation); %num of bits carried by one sample
num_bits = num_symbols*allocated_subchannels*bits_per_sample; % total number of bits transmitted
qam_sizes = [4 16 64 256];

cp_ratio = 0.25;%0.25;
cp_length = ceil(cp_ratio*M);

normalization = [2 10 42 170];

%--------Multiuser Parameters--------%
num_users = 5;
mapping = 1; % 0: LFDMA 1:DFDMA
start_indices = 10:2:28; %check if overlapping or overflowing
Q = ceil(M/N); % only for DFDMA mode bandwidth expansion factor

%--------Channel settings--------%
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
EbN0 = [20,10,2,15,10]; % SNR of the channel. noisy=1 to see the effects on channel
% SNR = [20,10,2,15]; % SNR of the channel. noisy=1 to see the effects on channel
noise_floor_dbHz = -165;
noise_floor_wattsHz = 10^(noise_floor_dbHz/10);
noise_pow = noise_floor_wattsHz*bw;

amplitudes = sqrt(noise_floor_wattsHz * 10.^(0.1*EbN0)*num_bits);
% amplitudes = [100,100,100,100,100]
user_sinr = zeros(1,num_users);
for i=1:num_users
    user_sinr(i) = 10*log10((amplitudes(i)^2)/(sum(amplitudes.^2)-(amplitudes(i)^2)+noise_pow));
end


%--------Parameter Check--------%
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

%%%%%%%%%%%%%%%% check boundaries



%--------Simulation Parameters--------%
EbN0_array = 0:1:100; % Array of SNR value(s) in dB.
save('EbN0_array.mat','EbN0_array');
BER=zeros(num_users,length(qam_sizes),length(EbN0_array));
ebno0_check_array = zeros(length(EbN0_array),length(qam_sizes));
num_trials = 1; % number of trials desired

if num_trials<1 || rem(num_trials,1)~=0
    num_trials
    error('Number of trials should be a positive integer.');
end

% if ~is_simulation
%     SNR = 1;
% end

ideal=0; %if channel is ideal

%% Additional configuration parameters
M_arr=M;
q_arr=qam_sizes;
s_arr=EbN0_array;
if is_simulation
    c1 = clock;
    conf=struct('M',M,...
        'mod_sch', q_arr,...
        'SNR_val',s_arr,...
        'num_users',num_users,...
        'user_indices',start_indices,...
        'mapping',mapping,...
        'num_symbols',num_symbols,...
        'num_trials',num_trials,...
        'started',c1,...
        'ended',[],...
        'mode','SC-FDMA',...
        'time_elapsed',0,...
        'noisy',noisy,...
        'explanation','SNR_val contains Eb/N0 values',...
        'version', 6);
    % save(sprintf('CONF%d-%d-%d-%d-%d.mat', c1(1:5)),'conf');
end

% 3- Print the configuration and ask for confirmation from user
disp('Configuration:')
if ~is_simulation
	disp(sprintf('N=%d\nM=%d\nnum_users=%d\nnum_symbols=%d\n%d-QAM\nnoisy=%d\nnum_bits=%d',N,M,num_users,num_symbols,modulation,noisy,num_bits))
	if noisy
	    disp(sprintf('Eb/N0=%d\n',EbN0(:)))
	    % disp(sprintf('fading=%d\n',fading(:)))
	    % disp(sprintf('profile=%s\n',profile{:}))
	else
	    disp(sprintf('SNR=N/A\n'))
	    % disp(sprintf('fading=%d\n',fading(:)))
	    % disp(sprintf('profile=%s\n',profile{:}))
	end
else
	conf
    disp(sprintf('Warning: All the BER/CONF/MSE files in current directory will be deleted!\n'))
    pause; 
    delete('BER*.mat');
    delete('CONF*.mat');
    conf.started = clock;
    save(sprintf('CONF%d-%d-%d-%d-%d.mat', conf.started(1:5)),'conf');
end