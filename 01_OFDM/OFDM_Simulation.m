%% OFDM_Simulation
%
% Burak Dayi
%
%

close all
clear all
clc
fprintf('--------------\n-----OFDM-----\n--simulation--\n--------------\n\n');
c1=clock;
fprintf('%d-%d-%d %d:%2d\n',c1(3),c1(2),mod(c1(1),100),c1(4),c1(5));
%% Transmission 
Config;
%disp('+Configuration is obtained.');
for M=M_arr
    for modulation=q_arr
        num_samples = M; %number of samples in a vector
        bits_per_sample = log2(modulation); %num of bits carried by one sample
        num_bits = num_symbols*num_samples*bits_per_sample; % total number of bits transmitted
        cp_length = ceil(cp_ratio*M);
        for SNR=s_arr
            if ~ideal
                disp(sprintf('M=%d, %d-PSK, SNR=%d dB, num_trials=%d, num_symbols=%d, num_bits=%d', M,modulation,SNR,num_trials,num_symbols,num_bits));
            else
                disp(sprintf('M=%d, %d-PSK, SNR=Ideal, num_trials=%d, num_symbols=%d, num_bits=%d', M,modulation,num_trials,num_symbols,num_bits));
            end
            for tt=1:num_trials
                Symbol_Creation;
                %disp('+Symbols are created.');
                Transmitter;
                %disp('+Transmitter Block is processed.');

                %% Channel
                Channel;
                %disp('+Channel effects are applied.');

                %% Reception
                Receiver;
                %disp('+Receiver Block is processed.');
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
    movefile(a(1,:),sprintf('CONF%d-%d-%d-%d-%d-Final.mat',c2(1:5)));
catch err
end

disp(sprintf('Simulation is completed in %f seconds', etime(c2,c1)));