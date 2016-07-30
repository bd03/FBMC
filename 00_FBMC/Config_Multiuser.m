%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Config_Multiuser.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created: 03-12-2014

% How to configure script?
%% 1- is_simulation setting: 
% set it 1 if you want to run a simulation with selected parameters. 
% set it 0 if you want to run only one try with selected parameters.
is_simulation = true;

%% 2- Set parameters for desired mode.
if is_simulation
    fprintf('--------------\n-----FBMC-----\n--simulation--\n--multiuser---\n--------------\n\n');
    mode_of_operation = 'SIMULATION';
    %% 2.a Simulation mode parameters
    %---- General filterbank parameters ----%
    K = 4; % overlapping factor    
    M = 256; % number of subchannels
    lp = K*M-1; % filter length
    delay = K*M+1-lp; %delay requirement
    symbol_delay = (K-1)+(K-1);
    
    % Prototype Filter coefficient generation
    h = func_Prototype_filter(M,K,lp);
    % disp('+Prototype filter has been designed.');

    %---- Prototype filter frequency coefficients----%
    %---------------------------------
    % This part should not be altered!
    %---------------------------------
    % K will select the row, the last column is background noise power in dB 
    % (that might be come in handy in future) for reference
    P=[ zeros(1,5);1 sqrt(2)/2 0 0 -35; 1 .911438 .411438 0 -44; 1 .97195983 sqrt(2)/2 .23514695 -65];
    % Prepare configuration containers:
    filter_params = struct;
    filter_params.K = K;
    filter_params.M = M;
    filter_params.lp = lp;
    filter_params.delay = delay;
    filter_params.symbol_delay = symbol_delay;
    filter_params.h = h;
    filter_params.P = P;

    %---- Supported values for simulation mode----%
    %---------------------------------
    % This part should not be altered!
    %---------------------------------
    % These are the values that are taken into account during the project.
    qam_sizes = [4 16 64 256]; % supported QAM modulation sizes
    EbN0_array = 0:1:40; % supported normalized EbN0 value(s) in dB.
    normalization = [2 10 42 170];
    % Prepare configuration containers:
    supported_vals = struct;
    supported_vals.qam_sizes = qam_sizes;
    supported_vals.EbN0_array = EbN0_array;
    supported_vals.normalization = normalization;


    %---- Simulation settings ----%
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
    autodefine = true;
    user_indices = [];
    if autodefine
        num_users = 5;
        subch_per_user = 20;
        % subch_per_user = M-2;
        user_indices = determine_indices(num_users,subch_per_user,M,2);
        allocated_subchannels = ones(1,num_users);
        allocated_subchannels = allocated_subchannels*subch_per_user;
    else
    %     user_indices = [1,50,170,220]; %case 1
        user_indices = [1,50,60,68,100,220,230,250]; %case 2
    %     user_indices = [1,50,60,68,170,220]; %case 3
    % user_indices = [20,200]; %case4
    user_indices=[1,150]
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
        num_users = length(allocated_subchannels);
        clear allocated_temp
    end

    num_frames = 10; % number of data frames in each FBMC block
    syms_per_frame = 20; %number of symbols per FBMC frame
    num_symbols = num_frames*syms_per_frame; % total number of data symbols
    num_trials = 1; % number of trials desired

    % vectors that hold information about number of bits and samples
    % num_symbols = num_frames*syms_per_frame; % total number of data symbols
    num_samples = allocated_subchannels*num_symbols; %number of samples in a vector
    % num_bits = num_samples*bits_per_sample; % total number of bits transmitted

    % M_arr=2.^(10:10); % array of M's that will be used in the simulation
    % q_arr=[4 16]; % array of QAM modes that will be used in sim.
    q_arr=[4 16 64 256]; % array of QAM modes that will be used in sim.
    % s_arr=200; % array of SNR values that will be used in the simulation
    s_arr=EbN0_array; % array of SNR values that will be used in the simulation
    % Prepare configuration containers
    sim_params = struct;
    sim_params.user_indices= user_indices;
    sim_params.allocated_subchannels= allocated_subchannels;
    sim_params.num_users= num_users;
    sim_params.num_frames= num_frames;
    sim_params.syms_per_frame= syms_per_frame;
    sim_params.num_symbols= num_symbols;
    sim_params.num_trials= num_trials;
    sim_params.num_samples= num_samples;
    sim_params.q_arr = q_arr;
    sim_params.s_arr = s_arr;
    
    %---- Channel settings ----%
    %rayleigh channel settings
    fading = 0; % set 0 for distortionless channel, set 1 for rayleigh channel
    % fading = 0; % set 0 for distortionless channel, set 1 for rayleigh channel
    use_matlab_channel = 0;
    bw = 5e6; % Transmission Bandwidth
    
    supported_channel_profiles = cell(3,1);
    supported_channel_profiles{1} = 'EPA';
    supported_channel_profiles{2} = 'EVA';
    supported_channel_profiles{3} = 'ETU'; % Valid channel profile selections
    ch_profile = NaN;
    if fading
        ch_profile = 'ETU';
    end

    ch_resp = NaN;
    if fading
        if use_matlab_channel
            max_doppler_shift = 1; %max. doppler shift
            [delay_a, pow_a] = LTE_channels2 (ch_profile,bw);
        %         ch_resp = rayleighchan(1/bw,max_doppler_shift,delay_a,pow_a); %channel model
            ch_resp = stdchan(1/bw,max_doppler_shift,'itur3GVAx'); % veh a
            if max_doppler_shift>0
                ch_resp.storeHistory = 1;
            end
            ch_resp.storePathGains =1;
        else
            ch_resp = LTE_channels (ch_profile,bw);        
        end
    end

    % noise settings
    % noisy = 0; %set 1 for EbN0 values to affect channel, set 0 for noiseless channel
    noisy = 1; %set 1 for EbN0 values to affect channel, set 0 for noiseless channel
    % EbN0 values for each user
    % EbN0 = [0,10]; %case 1
    % EbN0 = [10,10,10,10,10]; %case 2
    % EbN0 = [1,15,0]; %case 3
    % EbN0 = [10,20,15,-10]; %case2b
    % EbN0 = [4] %case4
    noise_floor_dbHz = -165;
    noise_floor_wattsHz = 10^(noise_floor_dbHz/10);
    noise_pow = noise_floor_wattsHz*bw;

    % % with the following addition, the system might be normalized in terms of bandwidth.
    % EbN0_array = EbN0_array + round(10*log10(bw));
    % Prepare configuration containers
    channel = struct;
    channel.fading = fading;
    channel.use_matlab_channel = use_matlab_channel;
    channel.bw = bw;
    channel.profile = ch_profile;
    channel.ch_resp = ch_resp;
    channel.noisy = noisy;
    channel.noise_floor_dbHz = noise_floor_dbHz;
    channel.noise_floor_wattsHz = noise_floor_wattsHz;
    channel.noise_pow = noise_pow;

    %---- Channel estimation settings ----%
    %---- Equalizer settings ----%
    eq_select = 2; % selection of equalizer type 1: one tap, 
    % 2: three tap w/ geometric interp, 3: three tap w/ linear interp
    % 4: no equalizer
    preamble_enabled = 1; %~(eq_select==4);

    %---- preamble setup ----%
    preamble = NaN;
    preamble_sel = NaN;
    extra_zero = NaN;
    zero_pads = NaN;
    fractional = NaN;
    % variables below automatically defined.
    length_preamble = NaN;
    estimation_column = NaN;

    if preamble_enabled
        % preamble_sel = 0: repmat([1 -j -1 j].',M/4,1) % IAM-I
        % preamble_sel = 1: repmat([1 1 -1 -1].',M/4,1) % IAM-R
        % preamble_sel = 2: repmat(repmat([1 -j -1 j].',M/4,1),1,3); % IAM-I with triple repetition.
        preamble_sel = 0;
        extra_zero = true; 
        zero_pads = 3; % 3ten az olunca bozuluyor---> extra zero eklenince 2'de de düzgün çıkıyor
        fractional = true; % if false, entire preamble will be applied, if true preamble samples at unused subchannels will be nulled.

        % everything below automatically defined
        [pr,ln,est]= func_preamble_creation(M, preamble_sel, zero_pads, extra_zero, user_indices, eq_select, fractional); 
        preamble = pr*sqrt(M/sum(allocated_subchannels));
        length_preamble = ln;
        estimation_column = est;
    end

    % Prepare configuration containers
    sp_params = struct;
    sp_params.eq_select = eq_select;
    sp_params.preamble_enabled = preamble_enabled;
    sp_params.preamble = preamble;
    sp_params.preamble_sel = preamble_sel;
    sp_params.extra_zero = extra_zero;
    sp_params.zero_pads = zero_pads;
    sp_params.fractional = fractional;
    sp_params.length_preamble = length_preamble;
    sp_params.estimation_column = estimation_column;

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
        q_arr
        error('Only 4-,16-,64-,128-,256-QAM schemes are supported. Define the number of constellation points (4,16,64,128,256) with modulation variable.')
    end
    
    % if ~all(ismember(M_arr, M_array)) || length(M_arr)<=0
    %     M_arr
    %     error('There are some values in M_array that are not supported.')
    % end
    
    if eq_select>4 || eq_select<1 || mod(eq_select,1)~=0
        eq_select
        error('eq_select should be an integer in range [1 4].');
    end
    
    % if ~ismember(valid_methods,estimation_method)
    %     estimation_method
    %     error('estimation_method should be either IAM, IAM4 or POP.');
    % end
    
    %---- Initialization of data containers ----%
    % BER matrix that will store BER values
    BER=zeros(num_users,length(qam_sizes),length(EbN0_array));
    % ebno0_check_array = zeros(length(EbN0_array),length(qam_sizes));
    % NMSE & MSE matrices that will store those values
    % MSE_f=zeros(length(M_array),length(qam_sizes),length(s_arr));
    % MSE_r=zeros(length(M_array),length(qam_sizes),length(s_arr));
    % MSE_a=zeros(length(M_array),length(qam_sizes),length(s_arr));

    % MSE_db_f=zeros(length(M_array),length(qam_sizes),length(s_arr));
    % MSE_db_r=zeros(length(M_array),length(qam_sizes),length(s_arr));
    % MSE_db_a=zeros(length(M_array),length(qam_sizes),length(s_arr));
    % CONF file that will store configuration of the simulation parameters
    c1 = clock; % time stamp
    conf=struct();
    % conf.M = M;
    conf=struct('filter_params', filter_params,...
        'sp_params', sp_params,...
        'supported_vals', supported_vals,...
        'sim_params',sim_params,...
        'channel',channel,...
        'started',c1,...
        'ended',[],...
        'time_elapsed',0,...
        'mode','FBMC',...
        'explanation','SNR_val contains Eb/N0 values',...
        'version', 9); % add later 'ch_profile', profile,...
    % conf.profile=profile;
    save(sprintf('CONF%d-%d-%d-%d-%d.mat', c1(1:5)),'conf');

    disp_config

    disp(sprintf('Press any key to proceed. \nIf you want to change configuration, please abort the script by pressing CTRL+C.'))
    disp(sprintf('Warning: All the BER/CONF/MSE files in current directory will be deleted!\n'))
    pause; 
    delete('BER*.mat');
    delete('CONF*.mat');
    conf.started = clock;
    save(sprintf('CONF%d-%d-%d-%d-%d.mat', conf.started(1:5)),'conf');
else
    fprintf('--------------\n-----FBMC-----\n--multiuser---\n--------------\n\n');
    mode_of_operation = 'MAIN';
    %% 2.b Main mode parameters
    %---- General filterbank parameters ----%
    K = 4; % overlapping factor    
    M = 256; % number of subchannels
    lp = K*M-1; % filter length
    delay = K*M+1-lp; %delay requirement
    symbol_delay = (K-1)+(K-1);

    %---- Prototype filter frequency coefficients----%
    %---------------------------------
    % This part should not be altered!
    %---------------------------------
    % K will select the row, the last column is background noise power in dB 
    % (that might be come in handy in future) for reference
    P=[ zeros(1,5);1 sqrt(2)/2 0 0 -35; 1 .911438 .411438 0 -44; 1 .97195983 sqrt(2)/2 .23514695 -65];
    % Prepare configuration containers:
    filter_params = struct;
    filter_params.K = K;
    filter_params.M = M;
    filter_params.lp = lp;
    filter_params.delay = delay;
    filter_params.symbol_delay = symbol_delay;
    filter_params.P = P;

    %---- Simulation settings ----%
    num_frames = 30; % number of data frames in each FBMC block
    syms_per_frame = 20; %number of symbols per FBMC frame
    modulation = 4; %4-, 16-, 64-, 256-QAM
    bits_per_sample = log2(modulation); %num of bits carried by one sample
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
    autodefine = true;
    user_indices = [];
    if autodefine
        num_users = 5;
        subch_per_user = 20;
        allocation_mode = 2;
        user_indices = determine_indices(num_users,subch_per_user,M,allocation_mode);
        allocated_subchannels = ones(1,num_users);
        allocated_subchannels = allocated_subchannels*subch_per_user;
    else
    %     user_indices = [1,50,170,220]; %case 1
        user_indices = [1,50,60,68,100,119,130,220,230,250]; %case 2
    %     user_indices = [1,50,60,68,170,220]; %case 3
    % user_indices = [1,M]; %case4
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
        num_users = length(allocated_subchannels);
        clear allocated_temp
    end

    % vectors that hold information about number of bits and samples
    num_symbols = num_frames*syms_per_frame; % total number of data symbols
    num_samples = allocated_subchannels*num_symbols; %number of samples in a vector
    num_bits = num_samples*bits_per_sample; % total number of bits transmitted
    % Prepare configuration containers
    sim_params = struct;
    sim_params.user_indices= user_indices;
    sim_params.allocated_subchannels= allocated_subchannels;
    sim_params.num_users= num_users;
    sim_params.num_frames= num_frames;
    sim_params.syms_per_frame= syms_per_frame;
    sim_params.num_symbols= num_symbols;
    sim_params.num_samples= num_samples;
    sim_params.num_bits = num_bits;
    sim_params.bits_per_sample = bits_per_sample;
    sim_params.modulation = modulation;
    % sim_params.amplitudes = defined after channel definition
    % sim_params.EbN0 = defined after channel definition

    %---- Channel settings ----%
    %rayleigh channel settings
    fading = repmat([1],num_users,1); % set 0 for distortionless channel, set 1 for rayleigh channel
    % fading = 0; % set 0 for distortionless channel, set 1 for rayleigh channel
    use_matlab_channel = 0;
    max_doppler_shift = 1;
    bw = 5e6; % Transmission Bandwidth
    
    supported_channel_profiles = cell(3,1);
    supported_channel_profiles{1} = 'EPA';
    supported_channel_profiles{2} = 'EVA';
    supported_channel_profiles{3} = 'ETU'; % Valid channel profile selections

    
    ch_profile = cell(num_users,1);
    for i=1:num_users
        ch_profile{i} = 'ETU';
        if ~fading(i)
            ch_profile{i} = NaN;
        end
    end

    ch_resp = NaN;
    if fading
        ch_resp = cell(num_users,1);

        % first create ch_resp for each profile
        % then fill in ch_resp cell.
        if use_matlab_channel
            % EPA
            [delay_a, pow_a] = LTE_channels2 ('EPA',bw);
            ch_resp_epa = rayleighchan(1/bw,max_doppler_shift,delay_a,pow_a); %channel model
            % ch_resp_epa = stdchan(1/bw,max_doppler_shift,'itur3GVAx'); % veh a
            if max_doppler_shift>0
                ch_resp_epa.storeHistory = 1;
            end
            ch_resp_epa.storePathGains =1;

            % EVA
            [delay_a, pow_a] = LTE_channels2 ('EVA',bw);
            ch_resp_eva = rayleighchan(1/bw,max_doppler_shift,delay_a,pow_a); %channel model
            % ch_resp_eva = stdchan(1/bw,max_doppler_shift,'itur3GVAx'); % veh a
            if max_doppler_shift>0
                ch_resp_eva.storeHistory = 1;
            end
            ch_resp_eva.storePathGains =1;

            % ETU
            [delay_a, pow_a] = LTE_channels2 ('ETU',bw);
            ch_resp_etu = rayleighchan(1/bw,max_doppler_shift,delay_a,pow_a); %channel model
            % ch_resp_etu = stdchan(1/bw,max_doppler_shift,'itur3GVAx'); % veh a
            if max_doppler_shift>0
                ch_resp_etu.storeHistory = 1;
            end
            ch_resp_etu.storePathGains =1;
        else
            % EPA
            ch_resp_epa = LTE_channels ('EPA',bw);
            % EVA
            ch_resp_eva = LTE_channels ('EVA',bw);
            % ETU
            ch_resp_etu = LTE_channels ('ETU',bw);
        end

        for i=1:num_users
            if ~isnan(ch_profile{i})
                if strcmp(ch_profile{i},'EPA')
                    ch_resp{i} = ch_resp_epa;
                elseif strcmp(ch_profile{i},'EVA')
                    ch_resp{i} = ch_resp_eva;
                elseif strcmp(ch_profile{i},'ETU')
                    ch_resp{i} = ch_resp_etu;
                else
                    error('Unknown channel profile')
                end
            else
                ch_resp{i} = NaN;
            end
        end
    end

    % noise settings
    % noisy = 0; %set 1 for EbN0 values to affect channel, set 0 for noiseless channel
    noisy = 1; %set 1 for EbN0 values to affect channel, set 0 for noiseless channel
    % EbN0 values for each user
    if noisy
        EbN0 = ones(1,num_users)*100; %case 2
        % EbN0 = [40, 50,60]
    else
        EbN0 = ones(1,num_users)*1000; % very high snr
    end        

    noise_floor_dbHz = -165;
    noise_floor_wattsHz = 10^(noise_floor_dbHz/10);
    noise_pow = noise_floor_wattsHz*bw;

    amplitudes = sqrt(noise_floor_wattsHz*(10.^(0.1*EbN0)).*num_bits);
    sim_params.amplitudes = amplitudes;
    sim_params.EbN0 = EbN0;
    % if noisy
    %     for i=1:num_users
    %         % required output power
    %         y_pow_req = noise_pow*10^(EbN0(i)/10);
    %         % current output power:
    %         % y_pow = mean(abs(y.*y)) %%normalized already  
    %         % adjust output power
    %         amplitudes(i) = sqrt(y_pow_req);
    %         % amplitudes(i) = sqrt(y_pow_req/y_pow)
    %     end        
    % else
    %     for i=1:num_users
    %         % required output power
    %         y_pow_req = noise_pow*10^(EbN0(i)/10);
    %         % current output power:
    %         % y_pow = mean(abs(y.*y)) %%normalized already  
    %         % adjust output power
    %         amplitudes(i) = 1;
    %     end        
    % end

    %plot allocation
    pl_y=zeros(1,M); % just in case
    color_array=['y', 'm', 'c', 'r', 'g', 'b', 'k'];
    pl_y(user_indices(1):user_indices(2)) = 1;
    stem(user_indices(1):user_indices(1), pl_y(user_indices(1):user_indices(1)),color_array(mod(1,length(color_array))+1));
    hold on
    for i=1:floor(length(user_indices)/2)
        pl_y(user_indices(2*i-1):user_indices(2*i)) = amplitudes(i);
        stem(user_indices(2*i-1):user_indices(2*i), pl_y(user_indices(2*i-1):user_indices(2*i)),color_array(mod(i,length(color_array))+1));
    end
    hold off
    axis([1 M 0 max(amplitudes)])
    xlabel('Subchannels')
    ylabel('Amplitude')
    title('Subchannel Allocation')

    % Prepare configuration containers
    channel = struct;
    channel.fading = fading;
    channel.use_matlab_channel = use_matlab_channel;
    channel.bw = bw;
    channel.profile = ch_profile;
    channel.ch_resp = ch_resp;
    channel.noisy = noisy;
    channel.noise_floor_dbHz = noise_floor_dbHz;
    channel.noise_floor_wattsHz = noise_floor_wattsHz;
    channel.noise_pow = noise_pow;


    %---- Channel estimation settings ----%
    %---- Equalizer settings ----%
    eq_select = 2; % selection of equalizer type 1: one tap, 
    % 2: three tap w/ geometric interp, 3: three tap w/ linear interp
    % 4: no equalizer
    preamble_enabled = ~(eq_select==4);

    %---- preamble setup ----%
    preamble = NaN;
    preamble_sel = NaN;
    extra_zero = NaN;
    zero_pads = NaN;
    fractional = NaN;
    % variables below automatically defined.
    length_preamble = NaN;
    estimation_column = NaN;

    if preamble_enabled
        % preamble_sel = 0: repmat([1 -j -1 j].',M/4,1) % IAM-I
        % preamble_sel = 1: repmat([1 1 -1 -1].',M/4,1) % IAM-R
        % preamble_sel = 2: repmat(repmat([1 -j -1 j].',M/4,1),1,3); % IAM-I with triple repetition.
        preamble_sel = 2;
        extra_zero = false;
        zero_pads = 1;
        fractional = true; % if false, entire preamble will be applied, if true preamble samples at unused subchannels will be nulled.

        % everything below automatically defined
        [pr,ln,est]= func_preamble_creation(M, preamble_sel, zero_pads, extra_zero, user_indices, eq_select, fractional); 
        preamble = pr; %*sqrt(M/sum(allocated_subchannels));
        length_preamble = ln;
        estimation_column = est;
    end

    % Prepare configuration containers
    sp_params = struct;
    sp_params.eq_select = eq_select;
    sp_params.preamble_enabled = preamble_enabled;
    sp_params.preamble = preamble;
    sp_params.preamble_sel = preamble_sel;
    sp_params.extra_zero = extra_zero;
    sp_params.zero_pads = zero_pads;
    sp_params.fractional = fractional;
    sp_params.length_preamble = length_preamble;
    sp_params.estimation_column = estimation_column;

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

    % if ~ismember(valid_methods,estimation_method)
    %     estimation_method
    %     error('estimation_method should be either IAM, IAM4 or POP.');
    % end

    if ~(length(EbN0)>=num_users)
        num_users
        EbN0
        error('Amplitude values should be defined for each user.');
    end

    % 3- Print the configuration and ask for confirmation from user
    disp_config

    disp(sprintf('Press any key to proceed. \nIf you want to change configuration, please abort the script by pressing CTRL+C.'))
    pause;

    disp('+Configuration is obtained.');
end