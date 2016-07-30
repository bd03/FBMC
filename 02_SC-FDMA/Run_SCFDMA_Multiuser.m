%% Run_SCFDMA_Multiuser
%
% Burak Dayi
%
% This is the main script that will run the multi-user case
% in SC-FDMA setup.
%
% Created: 05-02-2015

close all
clear all
clc

%% Configuration
Config_Multiuser;

if is_simulation
	fprintf('-----------------\n-----SC-FDMA-----\n---simulation----\n-----------------\n\n');
	c1=clock;
    fprintf('%d-%d-%d %d:%2d\n',c1(3),c1(2),mod(c1(1),100),c1(4),c1(5));

    % allocated_subchannels = N; %number of samples in a vector
    % preamble = [zeros(M,1) zeros(M,1) repmat([1 -j -1 j].',M/4,1) zeros(M,1) zeros(M,1)];
    % if strcmp(estimation_method,'IAM4')
    %     preamble =[zeros(M,1) preamble];
    % end 

    % conf.preamble=preamble;
    c1=clock;
    save(sprintf('CONF%d-%d-%d-%d-%d.mat', c1(1:5)),'conf');

    for modulation=q_arr
        bits_per_sample = log2(modulation); %num of bits carried by one sample
        num_bits = num_symbols*allocated_subchannels*bits_per_sample; % total number of bits transmitted
        for SNR=s_arr
        	amplitude = sqrt(noise_floor_wattsHz * 10^(0.1*SNR)*num_bits);
            if noisy
                disp(sprintf('M=%d, %d-QAM, Eb/N0=%d dB, num_trials=%d, num_symbols=%d, num_bits=%d', M,modulation,SNR,num_trials,num_symbols,num_bits));
            else
                disp(sprintf('M=%d, %d-QAM, Eb/N0=Ideal, num_trials=%d, num_symbols=%d, num_bits=%d', M,modulation,num_trials,num_symbols,num_bits));
            end
            for tt=1:num_trials
                % Symbol Creation
				symbols=cell(num_users,1);
				m=cell(num_users,1);
				bits=cell(num_users,1);
				for i=1:num_users
				    [m{i},bits{i},symbols{i}] = func_Symbol_Creation(num_bits,bits_per_sample, N, num_symbols,modulation);
				end
				% disp('+Symbols have been created.');

				% Transmitter
				tx_outputs=cell(num_users,1);
				A = zeros(num_users,1);
				out_without_ampl = cell(num_users,1);
				for i=1:num_users
				    [tx_outputs{i},A(i)] = func_Transmitter(N,M,Q,mapping,start_indices(i),symbols{i},cp_length,noisy,amplitude);
				end

				%% Channel
				ch_output = zeros(1,num_symbols*(M+cp_length));

				for i=1:num_users
					ch_output = ch_output + tx_outputs{i};
				end 

				if noisy
					noise = randn(1,num_symbols*(M+cp_length))*sqrt(noise_pow/2)+j*randn(1,num_symbols*(M+cp_length))*sqrt(noise_pow/2);
					ch_output = ch_output + noise;
				end

				mmm = find(qam_sizes==modulation);
			    snrrr = find(EbN0_array==SNR);
			    ebn0_check = 10*log10((sum(abs(tx_outputs{1}.*tx_outputs{1}))/num_bits)/noise_floor_wattsHz);
			    % ebn0_check = 10*log10((sum(abs(yy.*yy))/num_bits)/mean(abs(noise.*noise)));
			    ebno0_check_array(snrrr,mmm) = ebno0_check_array(snrrr,mmm)+(ebn0_check/num_trials);

				%!!! multipath channel effects could be implemented later.

				%% SC-FDMA Receiver
				% Receiver
				rx_output = cell(num_users,1);
				for i=1:num_users
					rx_output{i} = func_Receiver(N,M,Q,mapping,start_indices(i),ch_output,num_symbols,cp_length);
				end


				% Symbol Estimation
				bits_est=cell(num_users,1);
				m_est=cell(num_users,1);
				for i=1:num_users
					[bits_est{i},m_est{i}] = func_Symbol_Estimation(num_bits,bits_per_sample,N,rx_output{i},modulation,A(i));
				end

				%% Results
				% disp('+ Results are below');
				for i=1:num_users
				    % user = i
				    % if fading(i)
				    %     disp(sprintf('fading=%d profile=%s',fading(i),profile{i}))
				    % else
				    %     disp(sprintf('fading=%d profile=NaN',fading(i)))
				    % end
				    [serr,ser,berr,ber]=func_Results(m{i},m_est{i},bits{i},bits_est{i},num_symbols,num_bits,N);
				    BER(i,find(qam_sizes==modulation),find(EbN0_array==SNR)) = BER(i,find(qam_sizes==modulation),find(EbN0_array==SNR)) + ber/num_trials;
				    if noisy
				        % disp(sprintf('SNR=%d\n',SNR))
				    else
				    	% disp('SNR=N/A');
				    end

				    % disp(sprintf('Number of errorneous samples: %d/%d', serr,num_symbols*N))
					% disp(sprintf('SER: %f', ser));
					%     disp(sprintf('Number of errors qam_dif: %d', error2))
					%     disp(sprintf('Number of errors oqasdasdam_dif: %d', error3))
					% disp(sprintf('Number of bit errors: %d/%d', berr, num_bits))
					% disp(sprintf('BER: %f\n', ber))	    
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
            %print mean BER
            mean_ber = mean(BER(:,find(qam_sizes==modulation),find(EbN0_array==SNR)))
%             disp(mean_ber)
			if mean_ber == 0
                break
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
else
	% One-Run SC-FDMA Transmission
	%% SC-FDMA Transmitter
	% Symbol Creation
	symbols=cell(num_users,1);
	m=cell(num_users,1);
	bits=cell(num_users,1);
	for i=1:num_users
	    [m{i},bits{i},symbols{i}] = func_Symbol_Creation(num_bits,bits_per_sample, N, num_symbols,modulation);
	end
	disp('+Symbols have been created.');

	% Transmitter
	tx_outputs = cell(num_users,1);
	A = zeros(num_users,1);
	for i=1:num_users
	    [tx_outputs{i},A(i)] = func_Transmitter(N,M,Q,mapping,start_indices(i),symbols{i},cp_length,noisy,amplitudes(i));
	end

	%% Channel
	ch_output = zeros(1,num_symbols*(M+cp_length));

	for i=1:num_users
			ch_output = ch_output + tx_outputs{i};
		end 

	if noisy
		noise = randn(1,num_symbols*(M+cp_length))*sqrt(noise_pow/2)+j*randn(1,num_symbols*(M+cp_length))*sqrt(noise_pow/2);
		ch_output = ch_output + noise;
	end

	%!!! multipath channel effects could be implemented later.

	%% SC-FDMA Receiver
	% Receiver
	rx_output = cell(num_users,1);
	for i=1:num_users
		rx_output{i} = func_Receiver(N,M,Q,mapping,start_indices(i),ch_output,num_symbols,cp_length);
	end


	% Symbol Estimation
	bits_est=cell(num_users,1);
	m_est=cell(num_users,1);
	for i=1:num_users
	    [bits_est{i},m_est{i}] = func_Symbol_Estimation(num_bits,bits_per_sample,N,rx_output{i},modulation,A(i));
	end

	%% Results
	disp('+ Results are below');
	for i=1:num_users
	    user = i
	    % if fading(i)
	    %     disp(sprintf('fading=%d profile=%s',fading(i),profile{i}))
	    % else
	    %     disp(sprintf('fading=%d profile=NaN',fading(i)))
	    % end
	    [serr,ser,berr,ber]=func_Results(m{i},m_est{i},bits{i},bits_est{i},num_symbols,num_bits,N);
	    if noisy
	        disp(sprintf('Eb/N0=%d dB',EbN0(i)))
	        % disp(sprintf('Eb/=%f dB',user_sinr(i)))
	    else
	    	disp('SNR=N/A');
	    end

	    disp(sprintf('Number of errorneous samples: %d/%d', serr,num_symbols*N))
		disp(sprintf('SER: %f', ser));
		%     disp(sprintf('Number of errors qam_dif: %d', error2))
		%     disp(sprintf('Number of errors oqasdasdam_dif: %d', error3))
		disp(sprintf('Number of bit errors: %d/%d', berr, num_bits))
		disp(sprintf('BER: %f\n', ber))	    
	end
end