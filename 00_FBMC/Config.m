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
is_simulation = true;

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
    P=[ zeros(1,5);1 sqrt(2)/2 0 0 -35; 1 .911438 .411438 0 -44; 1 .97195983 sqrt(2)/2 .23514695 -65];
        
    %---- Supported values for simulation mode----%
    %---------------------------------
    % This part should not be altered!
    %---------------------------------
    % These are the values that are taken into account during the project.
    qam_sizes = [4 16 64 128 256]; % supported QAM modulation sizes
    SNR_array = 0:1:30; % supported SNR value(s) in dB.
    M_array = 2.^(2:10);
    normalization = [2 10 42 170];

    %---- Channel settings ----%
    % noise settings
    noisy = 1; %set 1 for SNR values to affect channel, set 0 for noiseless channel
    
    % rayleigh channel settings
    fading = 1; % set 0 for distortionless channel, set 1 for rayleigh channel
    bw = 5e+6; % Transmission Bandwidth
    max_doppler_shift = .0000001; %max. doppler shift in terms of hertz
    channel_profiles = ['EPA' 'EVA' 'ETU']; % Valid channel profile selections
    profile ='ETU'; %Channel profile
    use_matlab_channel = 0;
    if use_matlab_channel
%         [delay_a, pow_a] = LTE_channels2 (profile,bw);
%         ch_resp = rayleighchan(1/bw,max_doppler_shift,delay_a,pow_a); %channel model
        ch_resp = stdchan(1/bw,max_doppler_shift,'itur3GVAx'); % veh a
        if max_doppler_shift>0
            ch_resp.storeHistory = 1;
        end
        ch_resp.storePathGains =1;
    else
        ch_resp = LTE_channels (profile,bw);
    end
    
    %---- Channel estimation settings ----%
    % method
    valid_methods = ['IAM', 'IAM4']; % 'POP' may be added later.
    estimation_method = 'IAM';

    % preamble will be defined during the simulation    
    
    %---- Equalizer settings ----%
    eq_select = 2; % selection of equalizer type 1: one tap, 
    % 2: three tap w/ geometric interp, 3: three tap w/ linear interp
    % 4: no equalizer
    
    %---- Simulation settings ----%
    num_frames = 25; % number of data frames in each FBMC block
    syms_per_frame = 20; %number of symbols per FBMC frame
    num_symbols = num_frames*syms_per_frame; % total number of data symbols
    num_trials = 20; % number of trials desired
    
    M_arr=2.^(10:10); % array of M's that will be used in the simulation
    q_arr=[4 16 64 256]; % array of QAM modes that will be used in sim.
    s_arr=0:1:30; % array of SNR values that will be used in the simulation
        
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
    
    if ~ismember(valid_methods,estimation_method)
        estimation_method
        error('estimation_method should be either IAM, IAM4 or POP.');
    end
    
    %---- Initialization of data containers ----%
    % BER matrix that will store BER values
    BER=zeros(length(M_array),length(qam_sizes),length(s_arr));
    % NMSE & MSE matrices that will store those values
    MSE_f=zeros(length(M_array),length(qam_sizes),length(s_arr));
    MSE_r=zeros(length(M_array),length(qam_sizes),length(s_arr));
    MSE_a=zeros(length(M_array),length(qam_sizes),length(s_arr));

    MSE_db_f=zeros(length(M_array),length(qam_sizes),length(s_arr));
    MSE_db_r=zeros(length(M_array),length(qam_sizes),length(s_arr));
    MSE_db_a=zeros(length(M_array),length(qam_sizes),length(s_arr));
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
        'noisy',noisy,...
        'fading',fading,...
        'resp',ch_resp,...
        'ch_profile', profile,...
        'estimation_method',estimation_method,...
        'explanation','Blank',...
        'version', 5);
    save(sprintf('CONF%d-%d-%d-%d-%d.mat', c1(1:5)),'conf');
    
else
    %% 2.b Main mode parameters
    % 2.b Main mode parameters
    %---- General filterbank parameters ----%
    K = 4; % overlapping factor 
    M = 512; % number of subcarriers
    num_frames = 10; % number of data frames in each FBMC block
    syms_per_frame = 10; %number of symbols per FBMC frame
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
    P=[ zeros(1,5);1 sqrt(2)/2 0 0 -35; 1 .911438 .411438 0 -44; 1 .97195983 sqrt(2)/2 .23514695 -65];

    %---- Channel settings ----%
    % noise settings
    noisy = 0; %set 1 for SNR values to affect channel, set 0 for noiseless channel
    SNR = 10; % SNR of the channel. noisy=1 to see the effects on channel

    %rayleigh channel settings
    fading = 1; % set 0 for distortionless channel, set 1 for rayleigh channel
    bw = 5e+6; % Transmission Bandwidth
    max_doppler_shift = 1; %max. doppler shift
    channel_profiles = ['EPA', 'EVA', 'ETU']; % Valid channel profile selections
    profile ='EPA'; %Channel profile
    use_matlab_channel = 1;
    if use_matlab_channel
        [delay_a, pow_a] = LTE_channels2 (profile,bw);
%         ch_resp = rayleighchan(1/bw,max_doppler_shift,delay_a,pow_a); %channel model
        ch_resp = stdchan(1/bw,max_doppler_shift,'itur3GVAx'); % veh a
        if max_doppler_shift>0
            ch_resp.storeHistory = 1;
        end
        ch_resp.storePathGains =1;
    else
        ch_resp = LTE_channels (profile,bw);
    end

    %---- Channel estimation settings ----%
    % method
    valid_methods = ['IAM', 'IAM4']; % 'POP' may be added later.
    estimation_method = 'IAM';

    % preamble
    % IAM preambles 
    preamble = [zeros(M,1) repmat([1 1 -1 -1].',M/4,1) zeros(M,1)];
    % preamble = [zeros(M,1) repmat([1 -1 -1 1].',M/4,1) zeros(M,1)];
    % preamble = [zeros(M,1) repmat([1 -j -1 j].',M/4,1) zeros(M,1)];
%     preamble = [zeros(M,1) repmat([3 3 3*j 3*j -3*j -3*j -3 -3].',M/8,1) zeros(M,1)];
    % POP preambles
    % preamble = [repmat([1 -1].',M/2,1) zeros(M,1)];
    % IAM4 preambles
    % preamble = [zeros(M,1) zeros(M,1) repmat([1 1 -1 -1].',M/4,1) zeros(M,1)];
    if strcmp(estimation_method,'IAM4')
        preamble =[zeros(M,1) preamble];
    end

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

    if ~ismember(valid_methods,estimation_method)
        estimation_method
        error('estimation_method should be either IAM, IAM4 or POP.');
    end

end

%% 3- Print the configuration and ask for confirmation from user
disp('Configuration:')
if is_simulation
    conf
else
    if noisy
        disp(sprintf('K=%d, M=%d, num_symbols=%d, num_bits=%d, %d-QAM,\nnum_frames=%d, syms_per_frame=%d,fading=%d, equalizer=%d,\nestimation=%s, profile=%s, max_doppler_shift=%d, SNR=%d dB\nuse_matlab_channel=%d', K,M,num_symbols,num_bits,modulation,num_frames,syms_per_frame,fading,eq_select,estimation_method,profile,max_doppler_shift,SNR,use_matlab_channel));
    else
        disp(sprintf('K=%d, M=%d, num_symbols=%d, num_bits=%d, %d-QAM,\nnum_frames=%d, syms_per_frame=%d,fading=%d, equalizer=%d,\nestimation=%s, profile=%s, max_doppler_shift=%d, SNR=Ideal\nuse_matlab_channel=%d', K,M,num_symbols,num_bits,modulation,num_frames,syms_per_frame,fading,eq_select,estimation_method,profile,max_doppler_shift,use_matlab_channel));
    end
end

disp(sprintf('Press any key to proceed. \nIf you want to change configuration, please abort the script by pressing CTRL+C.'))

if is_simulation
    disp(sprintf('Warning: All the BER/CONF/MSE files in current directory will be deleted!\n'))
    pause; 
    delete('BER*.mat');
    delete('CONF*.mat');
    conf.started = clock;
    save(sprintf('CONF%d-%d-%d-%d-%d.mat', conf.started(1:5)),'conf');
else
    pause;
end











