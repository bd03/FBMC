%% zz_OneFileFBMC7.m
%
% Burak Dayi

% Experimental on Multiple Access Scheme


% Created: 21-11-2014

close all
clear all
clc
fprintf('--------------\n-----FBMC-----\n--------------\n\n');

%% Transmission 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Config
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 2.b Main mode parameters
%---- General filterbank parameters ----%
K = 4; % overlapping factor 
M = 256; % number of subcarriers
num_frames = 2000; % number of data frames in each FBMC block
syms_per_frame = 10; %number of symbols per FBMC frame
modulation = 4; %4-, 16-, 64-, 256-QAM
bits_per_sample = log2(modulation); %num of bits carried by one sample
lp = K*M-1; % filter length
% normalization = [2 10 42 170];
% q_arr=[4 16 64 256]; % array of QAM modes that will be used in sim.

%---- Multiple access settings----%
% Here the indices that are dedicated to different users can be
% defined.
% if autodefine is enabled, it would be sufficient to provide number of 
% users and number of subchannels allocated to each user. The index
% definitions are then automatically decided.
% In order to allocate different number of subchannels to different
% users, the indices should be manually defined. Note that for that
% sort of description, the autodefine variable must be disabled.
% !!!! the indices should be within [1,M]
autodefine = false;
user_indices = [];
if autodefine
    num_users = 5;
    subch_per_user = 30;
    user_indices = determine_indices(num_users,subch_per_user,M,1);
    allocated_subchannels = ones(1,num_users);
    allocated_subchannels = allocated_subchannels*subch_per_user;
else
    user_indices = [1,15,20,80,90,105,200,229];
    % quick parameter check
    % if the length is even
    if ~(mod(length(user_indices),2)==0)
        user_indices
        error('The length of user_indices vector should be even.');
    end
    
    % if the vector is monotonously increasing
    if ~(all(diff(user_indices)>0))
        user_indices
        error('The user_indices vector should be monotonously increasing.');
    end
    allocated_temp = diff(user_indices)+1;
    allocated_subchannels = allocated_temp(1:2:end);
    num_users = length(allocated_subchannels)
    clear allocated_temp
end

% plot allocation
pl_y=zeros(1,M); % just in case
color_array=['y', 'm', 'c', 'r', 'g', 'b', 'k'];
pl_y(user_indices(1):user_indices(2)) = 1;
stem(user_indices(1):user_indices(1), pl_y(user_indices(1):user_indices(1)),color_array(mod(1,length(color_array))+1));
hold on
for i=1:floor(length(user_indices)/2)
    pl_y(user_indices(2*i-1):user_indices(2*i)) = 1;
    stem(user_indices(2*i-1):user_indices(2*i), pl_y(user_indices(2*i-1):user_indices(2*i)),color_array(mod(i,length(color_array))+1));
end
hold off
axis([1 M 0 1.5])
xlabel('Subchannels')
title('Subchannel Allocation')

% vectors that hold information about number of bits and samples
num_symbols = num_frames*syms_per_frame; % total number of data symbols
num_samples = allocated_subchannels*num_symbols; %number of samples in a vector
num_bits = num_samples*bits_per_sample; % total number of bits transmitted

%---- Channel settings ----%
% noise settings
% noisy = 0; %set 1 for SNR values to affect channel, set 0 for noiseless channel
noisy = 1; %set 1 for SNR values to affect channel, set 0 for noiseless channel
SNR = 10; % SNR of the channel. noisy=1 to see the effects on channel
% SNR = 0; % SNR of the channel. noisy=1 to see the effects on channel

%rayleigh channel settings
fading = [1,1,0,0]; % set 0 for distortionless channel, set 1 for rayleigh channel
% fading = 0; % set 0 for distortionless channel, set 1 for rayleigh channel
bw = 5e+6; % Transmission Bandwidth
max_doppler_shift = 1; %max. doppler shift
channel_profiles = ['EPA', 'EVA', 'ETU']; % Valid channel profile selections
profile ='ETU'; %Channel profile
use_matlab_channel = 0;


%---- Channel estimation settings ----%
% method
valid_methods = ['IAM', 'IAM4']; % 'POP' may be added later.
estimation_method = 'IAM';

% preamble
% IAM preambles 
preamble = [zeros(M,1) repmat([1 -j -1 j].',M/4,1) zeros(M,1)];
% preamble = [zeros(M,1) zeros(M,1) repmat([1 -j -1 j].',M/4,1) zeros(M,1) zeros(M,1)];

% IAM4 preambles
% preamble = [zeros(M,1) zeros(M,1) repmat([1 1 -1 -1].',M/4,1) zeros(M,1)];
if strcmp(estimation_method,'IAM4')
    preamble =[zeros(M,1) preamble];
end

%---- Equalizer settings ----%
eq_select = 4; % selection of equalizer type 1: one tap, 
% 2: three tap w/ geometric interp, 3: three tap w/ linear interp
% 4: no equalizer

%---- Parameter check ---%
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

if ~~mod(log2(modulation),1)
    modulation
    error('Only modulation=2^m-QAM schemes are supported.');
end

if eq_select>4 || eq_select<1 || mod(eq_select,1)~=0
    eq_select
    error('eq_select should be an integer in range [1 4].');
end

if ~ismember(valid_methods,estimation_method)
    estimation_method
    error('estimation_method should be either IAM, IAM4 or POP.');
end


% 3- Print the configuration and ask for confirmation from user
disp('Configuration:')
if noisy
    disp(sprintf('K=%d, M=%d, num_symbols=%d, num_bits=%d, %d-QAM,\nnum_frames=%d, syms_per_frame=%d,fading=%d, equalizer=%d,\nestimation=%s, profile=%s, max_doppler_shift=%d, SNR=%d dB\nuse_matlab_channel=%d', K,M,num_symbols,num_bits,modulation,num_frames,syms_per_frame,fading,eq_select,estimation_method,profile,max_doppler_shift,SNR,use_matlab_channel));
else
    disp(sprintf('K=%d, M=%d, num_symbols=%d, num_bits=%d, %d-QAM,\nnum_frames=%d, syms_per_frame=%d,fading=%d, equalizer=%d,\nestimation=%s, profile=%s, max_doppler_shift=%d, SNR=Ideal\nuse_matlab_channel=%d', K,M,num_symbols,num_bits,modulation,num_frames,syms_per_frame,fading,eq_select,estimation_method,profile,max_doppler_shift,use_matlab_channel));
end
disp(sprintf('Press any key to proceed. \nIf you want to change configuration, please abort the script by pressing CTRL+C.'))
%pause;

disp('+Configuration is obtained.');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Prototype Filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h = func_Prototype_filter(M,K,lp);
disp('+Prototype filter is designed.');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Symbol_Creation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
symbols=cell(num_users,1);
m=cell(num_users,1);
bits=cell(num_users,1);
for i=1:num_users
    [m{i},bits{i},symbols{i}] = func_Symbol_Creation(num_bits(i),bits_per_sample, allocated_subchannels(i), num_symbols,modulation);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FBMC Transmitter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% OQAM Preprocessing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
oqam_m=cell(num_users,1);
for i=1:num_users
    [oqam_m{i},num_oqam_subsymbols] = func_OQAM_Preprocessing(M,num_frames,user_indices(2*i-1),user_indices(2*i),symbols{i},~(eq_select==4),preamble,syms_per_frame);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Synthesis Filter Bank
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
y=cell(num_users,1);
for i=1:num_users
    [y{i}] = func_Synthesis_Filter_Bank(M,K,h,num_oqam_subsymbols,oqam_m{i},lp);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Channel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
y_filtered=cell(num_users,1);
y_sum = zeros(1,K*M+(num_oqam_subsymbols-1)*M/2);
for i=1:num_users
%     y_ch{i} = func_Channel(y{i},fading,use_matlab_channel,profile,bw,noisy,SNR);
%     y_ch{i} = func_Channel(y{i},fading(i),use_matlab_channel,profile,bw,noisy(i),SNR(i));
    y_filtered{i} = func_Channel(y{i},fading(i),use_matlab_channel,profile,bw);
    y_sum = y_sum + y_filtered{i};
end

if ~noisy
    y_ch = y_sum;
else
    y_ch = awgn(y_sum,SNR,'measured'); %SNR in dB ~ 10log(Ps/Pn)
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FBMC Receiver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Analysis Filter Bank
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rx_output = func_Analysis_Filter_Bank(M,K,h,num_symbols,num_oqam_subsymbols,y_ch,lp);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Subchannel Processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sp_output = func_Subchannel_Processing(M,K,num_frames,syms_per_frame,rx_output,~(eq_select==4),preamble,estimation_method,eq_select);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% OQAM Postprocessing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
oqam_demod=cell(num_users,1);
for i=1:num_users
    oqam_demod{i} = func_OQAM_Postprocessing(M,num_symbols,user_indices(2*i-1),user_indices(2*i),sp_output);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Symbol Estmation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bits_est=cell(num_users,1);
m_est=cell(num_users,1);
for i=1:num_users
    [bits_est{i},m_est{i}] = func_Symbol_Estimation(allocated_subchannels(i),num_symbols,oqam_demod{i},modulation);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('+ Results are below');
for i=1:num_users
    user = i
    func_Results(m{i},m_est{i},bits{i},bits_est{i},num_symbols,num_bits(i),allocated_subchannels(i));
end


% close all

