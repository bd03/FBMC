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
    P=[ zeros(1,5);1 sqrt(2)/2 0 0 -35; 1 .911438 .411438 0 -44; 1 .971960 sqrt(2)/2 .235147 -65];
    
    %---- Supported values for simulation mode----%
    %---------------------------------
    % This part should not be altered!
    %---------------------------------
    % These are the values that are taken into account during the project.
    qam_sizes = [4 16 64 128 256]; % supported QAM modulation sizes
    SNR_array = 0:1:15; % supported SNR value(s) in dB.
    M_array = 2.^(2:9);
    
    %---- Simulation settings ----%
    num_symbols = 7; % number of symbols sent back to back in one transmission    
    num_trials = 7; % number of trials desired
    ideal=0; %set 0 for SNR values to affect channel, set 1 for ideal channel
    
    M_arr=2.^(4:6); % array of M's that will be used in the simulation
    q_arr=[4 16 64]; % array of QAM modes that will be used in sim.
    s_arr=SNR_array; % array of SNR values that will be used in the simulation
        
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
        'K',K,...
        'started',c1,...
        'ended',[],...
        'mode','FBMC',...
        'time_elapsed',0,...
        'ideal',ideal,...
        'resp',[],...
        'explanation','Blank');
    save(sprintf('CONF%d-%d-%d-%d-%d.mat', c1(1:5)),'conf');
    
else
    %% 2.b Main mode parameters
    %---- General filterbank parameters ----%
    K = 4; % overlapping factor
    M = 4; % number of subcarriers
    % num_frames = 0; % number of frames
    num_symbols = 7; % number of symbols sent back to back in one transmission
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
    
    %---- Channel settings ----%
    ideal=0; %set 0 for SNR values to affect channel, set 1 for ideal channel
    SNR = -10; % SNR of the channel. ideal=0 to see the effects on channel
    
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
        error('Only modulation=2^m-QAM schemes are supported.')
    end

end

%% 3- Print the configuration and ask for confirmation from user
disp('Configuration:')
if is_simulation
    conf
else
    if ~ideal
        disp(sprintf('K=%d, M=%d, num_symbols=%d, %d-QAM, num_bits=%d, SNR=%d dB', K,M,num_symbols,modulation,num_bits,SNR));
    else
        disp(sprintf('K=%d, M=%d, num_symbols=%d, %d-QAM, num_bits=%d, SNR=Ideal', K,M,num_symbols,modulation,num_bits));
    end
end
disp(sprintf('Press any key to proceed. \nIf you want to change configuration, please abort the script by pressing CTRL+C.'))
disp(sprintf('Warning: All the BER/CONF files in current directory will be deleted!\n'))
pause;
delete('BER*.mat');
delete('CONF*.mat');
conf.started = clock;
save(sprintf('CONF%d-%d-%d-%d-%d.mat', conf.started(1:5)),'conf');











