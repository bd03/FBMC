%% Run_FBMC
%
% Burak Dayi
%
% This is the main script that will replace FBMC_Simulation and FBMC_Main
% files. The tests can be run through this file.
%
% Created: 14-11-2014

close all
clear all
clc
%% Transmission 
Config;
if is_simulation
    fprintf('--------------\n-----FBMC-----\n--simulation--\n--------------\n\n');
    c1=clock;
    fprintf('%d-%d-%d %d:%2d\n',c1(3),c1(2),mod(c1(1),100),c1(4),c1(5));
    %paste FBMC_Simulation
    for M=M_arr
        lp = K*M-1; % filter length
        delay = K*M+1-lp; %delay requirement
        num_samples = M; %number of samples in a vector
        preamble = [zeros(M,1) zeros(M,1) repmat([1 -j -1 j].',M/4,1) zeros(M,1) zeros(M,1)];
        if strcmp(estimation_method,'IAM4')
            preamble =[zeros(M,1) preamble];
        end 

        conf.preamble=preamble;
        c1=clock;
        save(sprintf('CONF%d-%d-%d-%d-%d.mat', c1(1:5)),'conf');

        %Prototype filter
        Prototype_filter;

        for modulation=q_arr
            bits_per_sample = log2(modulation); %num of bits carried by one sample
            num_bits = num_symbols*num_samples*bits_per_sample; % total number of bits transmitted
            for SNR=s_arr
                if noisy
                    disp(sprintf('M=%d, %d-QAM, SNR=%d dB, num_trials=%d, num_symbols=%d, num_bits=%d', M,modulation,SNR,num_trials,num_symbols,num_bits));
                else
                    disp(sprintf('M=%d, %d-QAM, SNR=Ideal, num_trials=%d, num_symbols=%d, num_bits=%d', M,modulation,num_trials,num_symbols,num_bits));
                end
                for tt=1:num_trials
                    Symbol_Creation;
                    %disp('+Symbols are created.');
                    OQAM_Preprocessing;
                    %disp('+OQAM Preprocessing is done.');
                    Transmitter;
                    %disp('+Transmitter Block is processed.');

                    %% Channel
                    Channel;
                    %disp('+Channel effects are applied.');

                    %% Reception
                    Receiver;
                    %disp('+Receiver Block is processed.');
                    Subchannel_processing;
                    %disp('+Subchannel processing is done.');
                    OQAM_Postprocessing;
                    %disp('+OQAM Preprocessing is done.');
                    Symbol_Estimation;
                    %disp('+Symbol Estimation is done.');
                    Results;
                    %disp('+Calculations and Statistics..');
                end
            end
        end    
    end

    c2=clock;
    delete('BER*.mat');
    save(sprintf('BER%d-%d-%d-%d-%d-Final.mat', c2(1:5)),'BER');

    delete('MSE*.mat');
    try
        save(sprintf('MSE_f%d-%d-%d-%d-%d-Final.mat', c(1:5)),'MSE_f');
        save(sprintf('MSE_r%d-%d-%d-%d-%d-Final.mat', c(1:5)),'MSE_r');
        save(sprintf('MSE_a%d-%d-%d-%d-%d-Final.mat', c(1:5)),'MSE_a');

        save(sprintf('MSE_db_f%d-%d-%d-%d-%d-Final.mat', c(1:5)),'MSE_db_f');
        save(sprintf('MSE_db_r%d-%d-%d-%d-%d-Final.mat', c(1:5)),'MSE_db_r');
        save(sprintf('MSE_db_a%d-%d-%d-%d-%d-Final.mat', c(1:5)),'MSE_db_a');
    catch err
        warning('mse file');
    end


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
    %paste FBMC_Main
    fprintf('--------------\n-----FBMC-----\n--------------\n\n');
    disp('+Configuration is obtained.');
    Prototype_filter;
    disp('+Prototype filter is designed.');
    Symbol_Creation;
    disp('+Symbols are created.');
    OQAM_Preprocessing;
    disp('+OQAM Preprocessing is done.');
    Transmitter;
    disp('+Transmitter Block is processed.');

    %% Channel
    Channel;
    disp('+Channel effects are applied.');

    %% Reception
    Receiver;
    disp('+Receiver Block is processed.');
    Subchannel_processing;
    disp('+Subchannel processing is done.');
    OQAM_Postprocessing;
    disp('+OQAM Preprocessing is done.');
    Symbol_Estimation;
    disp('+Symbol Estimation is done.');
    Results;
    disp('+Calculations and Statistics..');
end