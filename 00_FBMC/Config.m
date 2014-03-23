%% Config
%
% Burak Dayi
%
% Configuraton file in which the user may set number of subcarriers, 
% overlapping factor, length of filter, prototype filter design parameters,
% modulation type, number of symbols (number of frames), SNR value and
% number of trials
%
% In the configuration file user input should be controlled and if
% necessary, error messages should be raised.
%
% Last updated: 18-03-2014
c1 = clock;

%% General filter parameters
K = 4; % overlapping factor
M = 4; % number of subcarriers
% num_frames = 0; % number of frames
num_symbols = 7; % number of symbols
num_samples = M; %number of samples in a vector
modulation = 4; %4-, 16-, 64-, 128-, 256-QAM
lp = K*M-1; % filter length
delay = K*M+1-lp; %delay requirement
bits_per_sample = log2(modulation); %num of bits carried by one sample
num_bits = num_symbols*num_samples*bits_per_sample; % total number of bits transmitted
qam_sizes = [4 16 64 128 256];
save('qam_sizes.mat','qam_sizes');



if K>4 || K<2
    K
    error('Range of integer K is [2 4]');
end

if ~~mod(log2(M),1)
    M
    error('The number of subcarriers M should be a power of 2.');
end

if ~(lp == K*M || lp == K*M-1 || lp == K*M+1)
    lp
    error('Only filter sizes of KM, KM+1, KM-1 are suported.');
end

if ~(ismember(modulation, qam_sizes))
    modulation
    %error('Only 4-QAM scheme is supported. Please set modulation = 4');
    error('Only 4-,16-,64-,128-,256-QAM schemes are supported. Define the number of constellation points (4,16,64,128,256) with modulation variable.')
end

%% Channel & simulation parameters
is_simulation = false;
SNR_array = 0:1:15; % Array of SNR value(s) in dB.
save('SNR_array.mat','SNR_array');
M_array = 2.^(2:9);
save('M_array.mat','M_array');
BER=zeros(length(M_array),length(qam_sizes),length(SNR_array));
num_trials = 5; % number of trials desired

if num_trials<1 || rem(num_trials,1)~=0
    num_trials
    error('Number of trials should be a positive integer.');
end

if ~is_simulation
    SNR = 10;
end

ideal=1; %if channel is ideal

%% Additional configuration parameters
M_arr=M_array;
q_arr=qam_sizes;
s_arr=1;
if is_simulation
    conf=struct('M_val',M_arr,'mod_sch', q_arr,'SNR_val',s_arr);
    save(sprintf('CONF%d-%d-%d-%d-%d.mat', c1(1:5)),'conf');
end

%% Prototype filter frequency coefficients
%---------------------------------
% This part should not be altered!
%---------------------------------
% K will select the row, the last column is background noise power in dB 
% (that might be come in handy in future) for reference
P=[ zeros(1,5);1 sqrt(2)/2 0 0 -35; 1 .911438 .411438 0 -44; 1 .971960 sqrt(2)/2 .235147 -65];
