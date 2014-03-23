%% Results
%
% Burak Dayi
%simulation
% This will present the results of simulation.
%
% Dependencies: m, m_est........
%
% Last updated: 18-03-2014

% disp('Results')

[error1,rate2] = symerr(m,m_est);
err=abs(m-m_est);
dif_qam =[qam_m qam_est];
error2 = symerr(qam_m,qam_est);
dif_oqam =[oqam_m oqam_demod];
error3 = symerr(oqam_m,oqam_demod);
dif_m = [m m_est];
[error4,rate] = symerr(bits,bits_est);

if is_simulation
    M_x= find(M_array==M);
    mod_x = find(qam_sizes==modulation);
    snr_x = find(SNR_array==SNR);
    BER(M_x, mod_x, snr_x) = BER(M_x, mod_x, snr_x)+rate/num_trials;
    if tt==num_trials
        disp(sprintf('BER = %f', BER(M_x, mod_x, snr_x)));
        c=clock;
        delete('BER*.mat');
        try
            save(sprintf('BER%d-%d-%d-%d-%d.mat', c(1:5)),'BER');
            
        catch err
            warning('BER of last iteration could not be saved.\n Do not worry, next iteration will already cover everthing.')
        end
        try
            a=ls('CONF2*');
            if ~strcmp(sprintf('CONF%d-%d-%d-%d-%d.mat',c(1:5)),a(1,:))
                movefile(a(1,:),sprintf('CONF%d-%d-%d-%d-%d.mat',c(1:5)));
            end
        catch err
            warning('conf file');
        end
    end
else
    if ~ideal
        disp(sprintf('\nK=%d, M=%d, num_sym=%d, %d-QAM, SNR=%d dB', K,M,num_symbols,modulation,SNR));
    else
        disp(sprintf('\nK=%d, M=%d, num_sym=%d, %d-QAM, SNR=Ideal', K,M,num_symbols,modulation));
    end
    disp(sprintf('Number of errorneous samples: %d/%d', error1,num_symbols*M))
    disp(sprintf('SER: %f', rate2));
%     disp(sprintf('Number of errors qam_dif: %d', error2))
%     disp(sprintf('Number of errors oqam_dif: %d', error3))
    disp(sprintf('Number of bit errors: %d/%d', error4,num_bits))
    disp(sprintf('BER: %f\n', rate))
end
