function [serr,ser,berr,ber] = func_Results(m,m_est,bits,bits_est,num_symbols,num_bits,allocated_subchannels)
%% func_Symbol_Creation
%
% Burak Dayi
%
% This function will present results for given data
%
% Created: 24-11-2014

%% Results

[serr,ser] = symerr(m,m_est);
err=abs(m-m_est);
bit_err = abs(bits-bits_est);
% dif_qam =[qam_m qam_est];
% error2 = symerr(qam_m(:,4:end),qam_est);
% dif_oqam =[oqam_m oqam_demod];
% error3 = symerr(oqam_m,oqam_demod);
% dif_m = [m m_est];
[berr,ber] = symerr(bits,bits_est);
% disp(sprintf('Number of errorneous samples: %d/%d', serr,num_symbols*allocated_subchannels))
% disp(sprintf('SER: %f', ser));
% %     disp(sprintf('Number of errors qam_dif: %d', error2))
% %     disp(sprintf('Number of errors oqasdasdam_dif: %d', error3))
% disp(sprintf('Number of bit errors: %d/%d', berr, num_bits))
% disp(sprintf('BER: %f\n', ber))
% 
% max_num_err=0;
% for i=1:num_symbols
%     num_errored(i) = length(find(err(:,i)));
%     num_errored_bits_per_symbol(i) = length(find(bit_err(1+(i-1)*num_bits/num_symbols:i*num_bits/num_symbols)));
%     
%     if max_num_err<num_errored(i)
%         max_num_err=num_errored(i);
%     end
% end    
% 
% errors = zeros(num_symbols,max_num_err);
% 
% for i=1:num_symbols
%     errors(i,1:length(find(err(:,i))))= find(err(:,i));
% end
% 
% for i=1:num_frames
%     num_errored_bits_per_frame(i) = length(find(bit_err(1+(i-1)*num_bits/num_frames:i*num_bits/num_frames)));
% end
% 
% if ~use_matlab_channel
%     % calculate MSE, if we use LTE channels
% %     pse_delta=zeros(1,M);
% %     pse_delta(1)=1;
% %     pse_delta_resp=conv(ch_resp,pse_delta);
%     fft_y=fft(ch_resp,M); % fft of impulse response
%         
%     %MSE
%     for i=1:num_frames
%         i
%         NMSE=(norm(fft_y.'-ch_resp_est(:,i)).^2)./(norm(fft_y.').^2);
%         NMSE_Db = 10*log10(NMSE)
%         
%         MSE=mean((abs(ch_resp_est(:,i)-fft_y.')).^2)
%         MSE_Db=10*log10(MSE)
%     end
%     
%     %plots
%     figure;
%     plot(abs(fft_y));
%     hold on
%     plot(abs(ch_resp_est(:,1)),'r');
%     % plot(abs(ch_resp_est(:,2)),'g.-');
%     % plot(abs(ch_resp_est(:,3)),'b.-');
%     
%     figure;
%     scatter(real(reshape(qam_est,num_symbols*M,1)),imag(reshape(qam_est,num_symbols*M,1)))
% else
%     figure;
%     plot(20*log10(abs(ch_resp_est(:,1))));
%     plot(ch_resp);
% end
% 
% 
% figfig=figure(54);
% subplot(311)
% plot(1:num_symbols,num_errored)
% xlabel('FBMC symbols');
% ylabel('Number of errors')
% title('Number of errors w.r.t symbols')
% subplot(312)
% plot(1:num_symbols,num_errored_bits_per_symbol)
% xlabel('FBMC symbols');
% ylabel('Number of errors')
% title('Number of bit errors w.r.t symbols')
% subplot(313)
% plot(1:num_frames,num_errored_bits_per_frame)
% xlabel('FBMC frames');
% ylabel('Number of errors')
% title('Number of bit errors w.r.t frames')
% 
% figu=figure(56);
% plot(1:num_symbols,errors,'o')
% xlabel('Symbols')
% ylabel('Subcarrier')
% title('Errors on each SC w.r.t symbols')
% c1=clock;
% savefig(figu,sprintf('%d%d%d%d%2d%d-figu',c1(3),c1(2),mod(c1(1),100),c1(4),c1(5),c1(6)*1000));
% savefig(figfig,sprintf('%d%d%d%d%2d%d-figfig',c1(3),c1(2),mod(c1(1),100),c1(4),c1(5),c1(6)*1000));
% save(sprintf('%d%d%d%d%2d%d-ch_resp.mat',c1(3),c1(2),mod(c1(1),100),c1(4),c1(5),c1(6)*1000),'ch_resp');
% save(sprintf('%d%d%d%d%2d%d-bits.mat',c1(3),c1(2),mod(c1(1),100),c1(4),c1(5),c1(6)*1000),'bits');