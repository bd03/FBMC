%% Run_FBMC_Multiuser
%
% Burak Dayi
%
% This is the main script that will replace FBMC_Simulation and FBMC_Main
% files. The tests can be run through this file.
%
% Created: 03-12-2014

close all
clear all
clc

% Configuration
Config_Multiuser;

if is_simulation
    c1=clock;
    fprintf('%d-%d-%d %d:%2d\n',c1(3),c1(2),mod(c1(1),100),c1(4),c1(5));
    %paste FBMC_Simulation
        
    c1=clock;
    save(sprintf('CONF%d-%d-%d-%d-%d.mat', c1(1:5)),'conf');

    for modulation=sim_params.q_arr
        bits_per_sample = log2(modulation); %num of bits carried by one sample
        num_bits = sim_params.num_samples*bits_per_sample; % total number of bits transmitted
        for s=sim_params.s_arr
            SNR = s + round(10*log10(bw));
            amplitudes = sqrt(channel.noise_floor_wattsHz*(10^(0.1*SNR)).*num_bits);
            % amplitudes = sqrt(channel.noise_floor_wattsHz *10*num_bits);
            if channel.noisy
                disp(sprintf('M=%d, %d-QAM, Eb/N0=%d dB, num_trials=%d, num_symbols=%d, total num_bits=%d', filter_params.M,modulation,s,sim_params.num_trials,sim_params.num_symbols,sum(num_bits)));
            else
                disp(sprintf('M=%d, %d-QAM, Eb/N0=Ideal, num_trials=%d, num_symbols=%d, total num_bits=%d', filter_params.M,modulation,sim_params.num_trials,sim_params.num_symbols,sum(num_bits)));
            end
            for tt=1:sim_params.num_trials
                %% FBMC Transmitter
                % Symbol_Creation
                symbols=cell(sim_params.num_users,1);
                m=cell(sim_params.num_users,1);
                bits=cell(sim_params.num_users,1);
                for i=1:sim_params.num_users
                    [m{i},bits{i},symbols{i}] = func_Symbol_Creation(num_bits(i),bits_per_sample, sim_params.allocated_subchannels(i), sim_params.num_symbols,modulation);
                end
                % disp('+Symbols have been created.');

                % OQAM Preprocessing
                oqam_m=cell(sim_params.num_users,1);
                for i=1:sim_params.num_users
                    [oqam_m{i},num_oqam_subsymbols] = func_OQAM_Preprocessing(symbols{i},filter_params.M,sim_params.num_frames,sim_params.user_indices(2*i-1),sim_params.user_indices(2*i),sp_params.preamble_enabled,sp_params.preamble,sim_params.syms_per_frame);
                end
                % disp('+OQAM Preprocessing has been completed.');

                % Synthesis Filter Bank

                y=cell(sim_params.num_users,1);
                A=zeros(sim_params.num_users,1); % constant value which is used to adjust amplitude of the output signal for each user
                for i=1:sim_params.num_users
                    [A(i),y{i}] = func_Synthesis_Filter_Bank(oqam_m{i},filter_params.M,filter_params.K,filter_params.h,num_oqam_subsymbols,filter_params.lp,channel.noisy,amplitudes(i));
                end
                % disp('+Synthesis Filter Bank has been processed.');

                %% Channel
                y_ch = func_Channel(y,channel,mode_of_operation);

                %% FBMC Receiver
                % Analysis Filter Bank
                rx_output = func_Analysis_Filter_Bank(y_ch,filter_params.M,filter_params.K,filter_params.h,sim_params.num_symbols,num_oqam_subsymbols,filter_params.lp,channel.noisy,A(i));

                % Subchannel Processing
                sp_output_cell=cell(sim_params.num_users,1); % keeps user data separately
                sp_output=zeros(filter_params.M,sim_params.num_symbols*2); % keeps all data together, matter of implementation
                for i=1:sim_params.num_users
                    [sp_output_cell{i}] = func_Subchannel_Processing(rx_output,filter_params.M,filter_params.K,sim_params.num_frames,sim_params.syms_per_frame,sim_params.user_indices(2*i-1),sim_params.user_indices(2*i),sp_params,sim_params.num_users);
                    % i
                    % size(sp_output_cell{i})
                    % sim_params.user_indices(2*i-1):sim_params.user_indices(2*i)
                    % sp_output_cell
                    % sim_params.user_indices(2*i-1):sim_params.user_indices(2*i)
                    % size(sp_output)
                    sp_output(sim_params.user_indices(2*i-1):sim_params.user_indices(2*i),:) = sp_output_cell{i};
                end

                % OQAM Postprocessing
                oqam_demod=cell(sim_params.num_users,1);
                for i=1:sim_params.num_users
                    oqam_demod{i} = func_OQAM_Postprocessing(sp_output,filter_params.M,sim_params.num_symbols,sim_params.user_indices(2*i-1),sim_params.user_indices(2*i));
                    % oqam_demod{i} = func_OQAM_Postprocessing(filter_params.M,sim_params.num_symbols,sim_params.user_indices(2*i-1),sim_params.user_indices(2*i),sp_output_cell{i});
                end

                % Symbol Estmation
                bits_est=cell(sim_params.num_users,1);
                m_est=cell(sim_params.num_users,1);
                for i=1:sim_params.num_users
                    [bits_est{i},m_est{i}] = func_Symbol_Estimation(oqam_demod{i},num_bits(i),sim_params.allocated_subchannels(i),sim_params.num_symbols,modulation);
                end

                % Results
                % disp('+ Results are below');
                for i=1:sim_params.num_users
                    % user = i
                    % if fading(i)
                    %     disp(sprintf('fading=%d profile=%s',fading(i),profile{i}))
                    % else
                    %     disp(sprintf('fading=%d profile=NaN',fading(i)))
                    % end
                    % if channel.noisy
                    %     disp(sprintf('Eb/N0=%d dB',EbN0(i)))
                    %     % disp(sprintf('SINR=%f dB', user_sinr(i)));
                    % else
                    %     disp('Eb/N0 = N/A')
                    % end
                    
                    [serr,ser,berr,ber]=func_Results(m{i},m_est{i},bits{i},bits_est{i},sim_params.num_symbols,num_bits(i),sim_params.allocated_subchannels(i));
                    BER(i,find(supported_vals.qam_sizes==modulation),find(supported_vals.EbN0_array==s)) = BER(i,find(supported_vals.qam_sizes==modulation),find(supported_vals.EbN0_array==s)) + ber/sim_params.num_trials;
                end
            end
            c=clock;
            delete('BER*.mat');
            try
                save(sprintf('BER%d-%d-%d-%d-%d.mat', c(1:5)),'BER');

            catch err
                warning('BER of last iteration could not be saved.\n Do not worry, next iteration will already cover everthing.')
            end
            
            try
                a=ls('CONF2*');
                if ~strcmp(version('-release'),'2013b')
                    a=a(1:end-1);
                end
                if ~strcmp(sprintf('CONF%d-%d-%d-%d-%d.mat',c(1:5)),a)
                    movefile(a,sprintf('CONF%d-%d-%d-%d-%d.mat',c(1:5)));
                end
            catch err
                warning('conf file');
            end
            mean_ber = mean(BER(:,find(supported_vals.qam_sizes==modulation),find(supported_vals.EbN0_array==s)))

            if mean_ber == 0
                break
            end

        end
    end    

    c2=clock;
    delete('BER*.mat');
    save(sprintf('BER%d-%d-%d-%d-%d-Final.mat', c2(1:5)),'BER');

    % delete('MSE*.mat');
    % try
    %     save(sprintf('MSE_f%d-%d-%d-%d-%d-Final.mat', c(1:5)),'MSE_f');
    %     save(sprintf('MSE_r%d-%d-%d-%d-%d-Final.mat', c(1:5)),'MSE_r');
    %     save(sprintf('MSE_a%d-%d-%d-%d-%d-Final.mat', c(1:5)),'MSE_a');

    %     save(sprintf('MSE_db_f%d-%d-%d-%d-%d-Final.mat', c(1:5)),'MSE_db_f');
    %     save(sprintf('MSE_db_r%d-%d-%d-%d-%d-Final.mat', c(1:5)),'MSE_db_r');
    %     save(sprintf('MSE_db_a%d-%d-%d-%d-%d-Final.mat', c(1:5)),'MSE_db_a');
    % catch err
    %     warning('mse file');
    % end

    try
        a=ls('CONF2*');
        if ~strcmp(version('-release'),'2013b')
            a=a(1:end-1);
        end
        load(a);
        conf.ended=c2;
        conf.time_elapsed = etime(conf.ended,conf.started);
        movefile(a,sprintf('CONF%d-%d-%d-%d-%d-Final.mat',conf.ended(1:5)));
        save(sprintf('CONF%d-%d-%d-%d-%d-Final.mat',conf.ended(1:5)),'conf');
    catch err
    end

    disp(sprintf('Simulation is completed in %f seconds', conf.time_elapsed));
    %showplot;
else
    %% FBMC Transmitter
    % Prototype Filter
    filter_params.h = func_Prototype_filter(filter_params.M,filter_params.K,filter_params.lp);
    disp('+Prototype filter has been designed.');

    % Symbol_Creation
    symbols=cell(sim_params.num_users,1);
    m=cell(sim_params.num_users,1);
    bits=cell(sim_params.num_users,1);
    for i=1:sim_params.num_users
        [m{i},bits{i},symbols{i}] = func_Symbol_Creation(sim_params.num_bits(i),sim_params.bits_per_sample, sim_params.allocated_subchannels(i), sim_params.num_symbols,sim_params.modulation);
    end
    disp('+Symbols have been created.');

    % OQAM Preprocessing
    oqam_m=cell(sim_params.num_users,1);
    for i=1:sim_params.num_users
        [oqam_m{i},num_oqam_subsymbols] = func_OQAM_Preprocessing(symbols{i},filter_params.M,sim_params.num_frames,sim_params.user_indices(2*i-1),sim_params.user_indices(2*i),sp_params.preamble_enabled,sp_params.preamble,sim_params.syms_per_frame);
    end
    disp('+OQAM Preprocessing has been completed.');

    % Synthesis Filter Bank

    y=cell(sim_params.num_users,1);
    % yy=cell(sim_params.num_users,1);
    A=zeros(sim_params.num_users,1); % constant value which is used to adjust amplitude of the output signal for each user
    for i=1:sim_params.num_users
        % [A(i),y{i}] = func_Synthesis_Filter_Bank(filter_params.M,filter_params.K,filter_params.h,num_oqam_subsymbols,oqam_m{i},filter_params.lp,channel.noisy,channel.noise_pow,sim_params.amplitudes(i));
        [A(i),y{i}] = func_Synthesis_Filter_Bank(oqam_m{i},filter_params.M,filter_params.K,filter_params.h,num_oqam_subsymbols,filter_params.lp,channel.noisy,sim_params.amplitudes(i));
    end
    disp('+Synthesis Filter Bank has been processed.');

    %% Channel
    % y_filtered=cell(sim_params.num_users,1);
    % y_sum = zeros(1,filter_params.K*filter_params.M+(num_oqam_subsymbols-1)*filter_params.M/2);
    % noise = randn(1,filter_params.K*filter_params.M+(num_oqam_subsymbols-1)*filter_params.M/2)*sqrt(noise_pow/2)+j*randn(1,filter_params.K*filter_params.M+(num_oqam_subsymbols-1)*filter_params.M/2)*sqrt(noise_pow/2);
    % for i=1:sim_params.num_users
    % %     y_ch{i} = func_Channel(y{i},fading,use_matlab_channel,profile,bw,channel.noisy,SNR);
    % %     y_ch{i} = func_Channel(y{i},fading(i),use_matlab_channel,profile,bw,channel.noisy(i),SNR(i));
    %     y_filtered{i} = func_Channel(y{i},fading(i),use_matlab_channel,profile{i},bw);
    %     y_sum = y_sum + y_filtered{i};
    % end 

    % if ~channel.noisy
    %     y_ch = y_sum;
    % else
    % %     y_ch = awgn(y_sum,SNR,'measured'); %SNR in dB ~ 10log(Ps/Pn)
    %     noise = randn(1,filter_params.K*filter_params.M+(num_oqam_subsymbols-1)*filter_params.M/2)*sqrt(noise_pow/2)+j*randn(1,filter_params.K*filter_params.M+(num_oqam_subsymbols-1)*filter_params.M/2)*sqrt(noise_pow/2);
    %     y_ch = y_sum + noise;
    % end
    %% Channel
    y_ch = func_Channel(y,channel,mode_of_operation);

     %% FBMC Receiver
    % Analysis Filter Bank
    rx_output = func_Analysis_Filter_Bank(y_ch,filter_params.M,filter_params.K,filter_params.h,sim_params.num_symbols,num_oqam_subsymbols,filter_params.lp,channel.noisy,A(i));

    % Subchannel Processing
    sp_output_cell=cell(sim_params.num_users,1); % keeps user data separately
    sp_output=zeros(filter_params.M,sim_params.num_symbols*2); % keeps all data together, matter of implementation
    for i=1:sim_params.num_users
        [sp_output_cell{i}] = func_Subchannel_Processing(rx_output,filter_params.M,filter_params.K,sim_params.num_frames,sim_params.syms_per_frame,sim_params.user_indices(2*i-1),sim_params.user_indices(2*i),sp_params,sim_params.num_users);
        % i
        % size(sp_output_cell{i})
        % sim_params.user_indices(2*i-1):sim_params.user_indices(2*i)
        % sp_output_cell
        % sim_params.user_indices(2*i-1):sim_params.user_indices(2*i)
        % size(sp_output)
        sp_output(sim_params.user_indices(2*i-1):sim_params.user_indices(2*i),:) = sp_output_cell{i};
    end

    % OQAM Postprocessing
    oqam_demod=cell(sim_params.num_users,1);
    for i=1:sim_params.num_users
        oqam_demod{i} = func_OQAM_Postprocessing(sp_output,filter_params.M,sim_params.num_symbols,sim_params.user_indices(2*i-1),sim_params.user_indices(2*i));
        % oqam_demod{i} = func_OQAM_Postprocessing(filter_params.M,sim_params.num_symbols,sim_params.user_indices(2*i-1),sim_params.user_indices(2*i),sp_output_cell{i});
    end

    % Symbol Estmation
    bits_est=cell(sim_params.num_users,1);
    m_est=cell(sim_params.num_users,1);
    for i=1:sim_params.num_users
        [bits_est{i},m_est{i}] = func_Symbol_Estimation(oqam_demod{i},sim_params.num_bits(i),sim_params.allocated_subchannels(i),sim_params.num_symbols,sim_params.modulation);
    end

    % Results
    disp('+ Results are below');
    for i=1:sim_params.num_users
        user = i
        if channel.fading(i)
            disp(sprintf('fading=%d profile=%s',channel.fading(i),channel.profile{i}))
        else
            disp(sprintf('fading=%d profile=NaN',channel.fading(i)))
        end
        if channel.noisy
            disp(sprintf('Eb/N0=%d dB',sim_params.EbN0(i)))
            % disp(sprintf('SINR=%f dB', user_sinr(i)));
        else
            disp('Eb/N0 = N/A')
        end
        
        [serr,ser,berr,ber]=func_Results(m{i},m_est{i},bits{i},bits_est{i});

        disp(sprintf('Number of errorneous samples: %d/%d', serr,sim_params.num_symbols*sim_params.allocated_subchannels(i)))
        disp(sprintf('SER: %f', ser));
        %     disp(sprintf('Number of errors qam_dif: %d', error2))
        %     disp(sprintf('Number of errors oqasdasdam_dif: %d', error3))
        disp(sprintf('Number of bit errors: %d/%d', berr, sim_params.num_bits(i)))
        disp(sprintf('BER: %f\n', ber))    
    end
end