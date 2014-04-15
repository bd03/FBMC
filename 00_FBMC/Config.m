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
% Created: 02-03-2014

% How to configure script?
%% 1- is_simulation setting: 
% set it 1 if you want to run a simulation with selected parameters. 
% set it 0 if you want to run only one try with selected parameters.
is_simulation = false;

%% 2- Set parameters for desired mode.
if is_simulation
    %% 2.a Simulation mode parameters
    %---- General filterbank parameters ----%
    K = 4; % overlapping factor    
    %---- Prototype filter frequency coefficients----%
    %---------------------------------
    % This part should not be altered!
    %---------------------------------
    % K will select the row, the last column is background noise power in dB 
    % (that might be come in handy in future) for reference
    P=[ zeros(1,5);1 sqrt(2)/2 0 0 -35; 1 .911438 .411438 0 -44; 1 .971960 sqrt(2)/2 .235147 -65];
        
    %---- Supported values for simulation mode----%
    %---------------------------------
    % This part should not be altered!
    %---------------------------------
    % These are the values that are taken into account during the project.
    qam_sizes = [4 16 64 128 256]; % supported QAM modulation sizes
    SNR_array = 0:1:15; % supported SNR value(s) in dB.
    M_array = 2.^(2:10);

    %---- Channel settings ----%
    % noise settings
    noisy = 0; %set 1 for SNR values to affect channel, set 0 for noiseless channel
    
    % rayleigh channel settings
    fading = 1; % set 0 for distortionless channel, set 1 for rayleigh channel
    bw = 5e+6; % Transmission Bandwidth
    max_doppler_shift = 10; %max. doppler shift in terms of hertz
    channel_profiles = ['EPA' 'EVA' 'ETU']; % Valid channel profile selections
    profile ='EPA'; %Channel profile
    [delay_a pow_a] = LTE_channels (profile,bw);
    ch_resp = rayleighchan(1/bw,max_doppler_shift,delay_a,pow_a); %channel model
    ch_resp.storeHistory = 1;
    ch_resp.storePathGains =1;

    %---- Equalizer settings ----%
    eq_select = 2; % selection of equalizer type 1: one tap, 
    % 2: three tap w/ geometric interp, 3: three tap w/ linear interp
    % 4: no equalizer
    
    %---- Simulation settings ----%
    num_frames = 20; % number of data frames in each FBMC block
    syms_per_frame = 10; %number of symbols per FBMC frame
    num_symbols = num_frames*syms_per_frame; % total number of data symbols
    num_trials = 3; % number of trials desired
    
    M_arr=2.^(9); % array of M's that will be used in the simulation
    q_arr=[4]; % array of QAM modes that will be used in sim.
    s_arr=[1]; % array of SNR values that will be used in the simulation
        
    %---- Parameter check ----%
    if K>4 || K<2
        K
        error('Range of integer K is [2 4]');
    end
    
    if num_trials<1 || rem(num_trials,1)~=0
        num_trials
        error('Number of trials should be a positive integer.');
    end
    
    if ~all(ismember(q_arr, qam_sizes)) || length(q_arr)<=0
        modulation
        error('Only 4-,16-,64-,128-,256-QAM schemes are supported. Define the number of constellation points (4,16,64,128,256) with modulation variable.')
    end
    
    if ~all(ismember(M_arr, M_array)) || length(M_arr)<=0
        M_arr
        error('There are some values in M_array that are not supported.')
    end
    
    if eq_select>4 || eq_select<1 || mod(eq_select,1)~=0
        eq_select
        error('eq_select should be an integer in range [1 4].');
    end
    
    %---- Initialization of data containers ----%
    % BER matrix that will store BER values
    BER=zeros(length(M_array),length(qam_sizes),length(SNR_array));
    % CONF file that will store configuration of the simulation parameters
    c1 = clock; % time stamp
    conf=struct('M_val',M_arr,...
        'mod_sch', q_arr,...
        'SNR_val',s_arr,...
        'num_symbols',num_symbols,...
        'num_trials',num_trials,...
        'num_frames',num_frames,...
        'syms_per_frame',syms_per_frame,...
        'eq_select', eq_select,...
        'K',K,...
        'started',c1,...
        'ended',[],...
        'mode','FBMC',...
        'time_elapsed',0,...
        'ideal',noisy,...
        'resp',ch_resp,...
        'ch_profile', profile,...
        'explanation','Blank',...
        'version', 3);
    save(sprintf('CONF%d-%d-%d-%d-%d.mat', c1(1:5)),'conf');
    
else
    %% 2.b Main mode parameters
    % 2.b Main mode parameters
    %---- General filterbank parameters ----%
    K = 4; % overlapping factor 
    M = 512; % number of subcarriers
    num_frames = 20; % number of data frames in each FBMC block
    syms_per_frame = 30; %number of symbols per FBMC frame
    num_symbols = num_frames*syms_per_frame; % total number of data symbols
    num_samples = M; %number of samples in a vector
    modulation = 4; %4-, 16-, 64-, 128-, 256-QAM
    bits_per_sample = log2(modulation); %num of bits carried by one sample
    num_bits = num_symbols*num_samples*bits_per_sample; % total number of bits transmitted
    lp = K*M-1; % filter length
    delay = K*M+1-lp; %delay requirement

    %---- Prototype filter frequency coefficients----%
    %---------------------------------
    % This part should not be altered!
    %---------------------------------
    % K will select the row, the last column is background noise power in dB 
    % (that might be come in handy in future) for reference
    P=[ zeros(1,5);1 sqrt(2)/2 0 0 -35; 1 .911438 .411438 0 -44; 1 .971960 sqrt(2)/2 .235147 -65];

    %---- Preamble creation ----%
    preamble = [repmat([1 1 -1 -1].',M/4,1) zeros(M,1)]; %IAM
    %preamble = [repmat([1 -1].',M/2,1) zeros(M,1)]; %POP
    %preamble = [repmat([1 -j -1 j].',M/4,1) zeros(M,1)];
    %preamble = [repmat([1 -1 -1 1].',M/4,1) zeros(M,1)]; %IAM
    %preamble =  [repmat([-3 -3 -1 -1 1 1 3 3].',M/8,1) zeros(M,1)];

    %---- Channel settings ----%
    % noise settings
    noisy = 0; %set 1 for SNR values to affect channel, set 0 for noiseless channel
    SNR = 20; % SNR of the channel. ideal=0 to see the effects on channel

    %rayleigh channel settings
    fading = 1; % set 0 for distortionless channel, set 1 for rayleigh channel
    bw = 5e+6; % Transmission Bandwidth
    max_doppler_shift = 10; %max. doppler shift in terms of hertz
    channel_profiles = ['EPA' 'EVA' 'ETU']; % Valid channel profile selections
    profile ='EPA'; %Channel profile
    [delay_a pow_a] = LTE_channels (profile,bw);
    ch_resp = rayleighchan(1/bw,max_doppler_shift,delay_a,pow_a); %channel model
    ch_resp.storeHistory = 1;
    ch_resp.storePathGains =1;

    %---- Equalizer settings ----%
    eq_select = 2; % selection of equalizer type 1: one tap, 
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

end

%% 3- Print the configuration and ask for confirmation from user
disp('Configuration:')
if is_simulation
    conf
else
    if noisy
        disp(sprintf('K=%d, M=%d, num_symbols=%d, %d-QAM, num_bits=%d,\nnum_frames=%d, syms_per_frames=%d,\nfading=%d, equalizer=%d, estimation=''POP'', SNR=%d dB', K,M,num_symbols,modulation,num_bits,num_frames,syms_per_frame,fading,eq_select,SNR));
    else
        disp(sprintf('K=%d, M=%d, num_symbols=%d, %d-QAM, num_bits=%d,\nnum_frames=%d, syms_per_frames=%d,\nfading=%d, equalizer=%d, estimation=''POP'', SNR=Ideal', K,M,num_symbols,modulation,num_bits,num_frames,syms_per_frame,fading,eq_select));
    end
end

disp(sprintf('Press any key to proceed. \nIf you want to change configuration, please abort the script by pressing CTRL+C.'))
disp(sprintf('Warning: All the BER/CONF files in current directory will be deleted!\n'))
pause; 
delete('BER*.mat');
delete('CONF*.mat');
conf.started = clock;
save(sprintf('CONF%d-%d-%d-%d-%d.mat', conf.started(1:5)),'conf');











