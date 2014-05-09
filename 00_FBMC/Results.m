%% Results
%
% Burak Dayi
%simulation
% This will present the results of simulation.
%
% Dependencies: m, m_est........
%
% Created: 02-03-2014

% disp('Results')

[error1,rate2] = symerr(m,m_est);
err=abs(m-m_est);
bit_err = abs(bits-bits_est);
% dif_qam =[qam_m qam_est];
% error2 = symerr(qam_m,qam_est);
% dif_oqam =[oqam_m oqam_demod];
% error3 = symerr(oqam_m,oqam_demod);
% dif_m = [m m_est];
[error4,rate] = symerr(bits,bits_est);

if is_simulation
    M_x= find(M_array==M);
    mod_x = find(qam_sizes==modulation);
    snr_x = find(s_arr==SNR);
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
            if ~strcmp(version('-release'),'2013b')
                a=a(1:end-1);
            end
            if ~strcmp(sprintf('CONF%d-%d-%d-%d-%d.mat',c(1:5)),a)
                movefile(a,sprintf('CONF%d-%d-%d-%d-%d.mat',c(1:5)));
            end
        catch err
            warning('conf file');
        end
    end
else
    
    disp(sprintf('Number of errorneous samples: %d/%d', error1,num_symbols*M))
    disp(sprintf('SER: %f', rate2));
    %     disp(sprintf('Number of errors qam_dif: %d', error2))
    %     disp(sprintf('Number of errors oqasdasdam_dif: %d', error3))
    disp(sprintf('Number of bit errors: %d/%d', error4,num_bits))
    disp(sprintf('BER: %f\n', rate))

    max_num_err=0;
    for i=1:num_symbols
        num_errored(i) = length(find(err(:,i)));
        num_errored_bits_per_symbol(i) = length(find(bit_err(1+(i-1)*num_bits/num_symbols:i*num_bits/num_symbols)));

        if max_num_err<num_errored(i)
            max_num_err=num_errored(i);
        end
    end    

    errors = zeros(num_symbols,max_num_err);

    for i=1:num_symbols
        errors(i,1:length(find(err(:,i))))= find(err(:,i));
    end

    for i=1:num_frames
        num_errored_bits_per_frame(i) = length(find(bit_err(1+(i-1)*num_bits/num_frames:i*num_bits/num_frames)));
    end

    figfig=figure(54);
    subplot(311)
    plot(1:num_symbols,num_errored)
    xlabel('FBMC symbols');
    ylabel('Number of errors')
    title('Number of errors w.r.t symbols')
    subplot(312)
    plot(1:num_symbols,num_errored_bits_per_symbol)
    xlabel('FBMC symbols');
    ylabel('Number of errors')
    title('Number of bit errors w.r.t symbols')
    subplot(313)
    plot(1:num_frames,num_errored_bits_per_frame)
    xlabel('FBMC frames');
    ylabel('Number of errors')
    title('Number of bit errors w.r.t frames')

    figu=figure(56);
    plot(1:num_symbols,errors,'o')
    xlabel('Symbols')
    ylabel('Subcarrier')
    title('Errors on each SC w.r.t symbols')
    c1=clock;
    % savefig(figu,sprintf('%d%d%d%d%2d%d-figu',c1(3),c1(2),mod(c1(1),100),c1(4),c1(5),c1(6)*1000));
    % savefig(figfig,sprintf('%d%d%d%d%2d%d-figfig',c1(3),c1(2),mod(c1(1),100),c1(4),c1(5),c1(6)*1000));
    % save(sprintf('%d%d%d%d%2d%d-ch_resp.mat',c1(3),c1(2),mod(c1(1),100),c1(4),c1(5),c1(6)*1000),'ch_resp');
    % save(sprintf('%d%d%d%d%2d%d-bits.mat',c1(3),c1(2),mod(c1(1),100),c1(4),c1(5),c1(6)*1000),'bits');
end
