%% FBMC_Simulation
%
% Burak Dayi
%
% This is the simulation script that will call other scripts/functions.
%
% This will implement Filterbank Multicarrier transmission and reception
% with modulation and demodulations steps described in PHYDYAS documents
% Deliverable 5.1 and FBMC Physical Layer: A primer.
%
% Notes:
% Configuration will allow only subcarrier sizes of a power of two so as to
% be able to use FFT in the implementation
%
% The prototype filter used here is the PHYDYAS filter described in
% Deliverable 5.1: Prototype filter and structure optimization
% The implementation of the blocks follows what is described in Deliverable
% 5.1 transmultiplexer architecture.
%
% Channel implementation
% Subchannel processing
% Calculations and statistics
%
%!!!!!!!!!!!!!!!!!!clear all configdeki deðiþiklikler
% Last updated: 17-03-2014 

close all
% clear all
clc
fprintf('--------------\n-----FBMC-----\n--simulation--\n--------------\n\n');
c1=clock;
fprintf('%d-%d-%d %d:%2d\n',c1(3),c1(2),mod(c1(1),100),c1(4),c1(5));
%% Transmission 
Config;
%disp('+Configuration is obtained.');
for M=M_arr
    lp = K*M-1; % filter length
    delay = K*M+1-lp; %delay requirement
    num_samples = M; %number of samples in a vector
    
    Prototype_filter;
    %disp('+Prototype filter is designed.');
    for modulation=q_arr
        bits_per_sample = log2(modulation); %num of bits carried by one sample
        num_bits = num_symbols*num_samples*bits_per_sample; % total number of bits transmitted
        for SNR=s_arr
            if ~ideal
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
try
    a=ls('CONF2*'); 
    load(a);
    conf.ended=c2;
    conf.time_elapsed = etime(conf.ended,conf.started);
    movefile(a(1,:),sprintf('CONF%d-%d-%d-%d-%d-Final.mat',conf.ended(1:5)));
    save(sprintf('CONF%d-%d-%d-%d-%d-Final.mat',conf.ended(1:5)),'conf');
catch err
end

disp(sprintf('Simulation is completed in %f seconds', conf.time_elapsed));
%showplot;