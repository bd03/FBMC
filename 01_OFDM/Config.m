%% Config
%
% Burak Dayi
%
% Configuraton file in which the user may set number of subcarriers, 
% modulation type, number of symbols (number of frames), SNR value and
% number of trials
%
% In the configuration file user input should be controlled and if
% necessary, error messages should be raised.
%
% Last updated: 18-03-2014
c1 = clock;

%% General parameters
M = 128; % number of subcarriers
num_symbols = 5000; % number of symbols
modulation = 4; %4-, 16-, 64-, 128-, 256-QAM
num_samples = M; %number of samples in a vector
bits_per_sample = log2(modulation); %num of bits carried by one sample
num_bits = num_symbols*num_samples*bits_per_sample; % total number of bits transmitted
qam_sizes = [2 4];

cp_ratio = 0.25;
cp_length = ceil(cp_ratio*M);

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
is_simulation = true;
SNR_array = 0:1:15; % Array of SNR value(s) in dB.
save('SNR_array.mat','SNR_array');
M_array = 2.^(2:9);
save('M_array.mat','M_array');
BER=zeros(length(M_array),length(qam_sizes),length(SNR_array));
num_trials = 200; % number of trials desired

if num_trials<1 || rem(num_trials,1)~=0
    num_trials
    error('Number of trials should be a positive integer.');
end

if ~is_simulation
    SNR = 1;
end

ideal=0; %if channel is ideal

%% Additional configuration parameters
M_arr=M_array;
q_arr=qam_sizes;
s_arr=SNR_array;
if is_simulation
    conf=struct('M_val',M_arr,'mod_sch', q_arr,'SNR_val',s_arr);
    save(sprintf('CONF%d-%d-%d-%d-%d.mat', c1(1:5)),'conf');
end