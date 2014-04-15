%% Subchannel_processing
%
% Burak Dayi
%
% This script will perform subchannel processing steps.
%
% Dependencies: rx_output - receiver block output signal
% Output: sp_output - processed signal output
%
% Created: 02-03-2014

% disp('Subchannel Processing')

% remove preamble
scp_preamble = rx_output(:,2*K-1:2*(syms_per_frame+1):end-(2*K-1));
scp_data = rx_output;

% %%!!!!!!!!!!!!!!!!hack!!!!!!!!!!!!!!!!!!!!!!!!!!
num_symbols = num_symbols-num_frames; %num_frames= no of preamble
% 
% %estimation of channel
for i=1:num_frames
    ch_resp_est(:,i) = scp_preamble(:,i)./(preamble(:,1)*sumfactor);
end

if eq_select == 1
    % one tap equalizer
    for ii=1:num_frames
        sp_output(:,1+(ii-1)*2*syms_per_frame:2*syms_per_frame+(ii-1)*2*syms_per_frame)...
            = scp_data(:,2*K+1+(ii-1)*2*syms_per_frame+2*(ii-1):2*K+2*syms_per_frame+(ii-1)*2*syms_per_frame+2*(ii-1));
        for iii=1:2*syms_per_frame
            sp_output(:,1+(ii-1)*2*syms_per_frame+iii-1)=sp_output(:,1+(ii-1)*2*syms_per_frame+iii-1)./ch_resp_est(:,ii);
        end
    end
elseif eq_select == 2 || eq_select == 3
    % % three tap equalizer
    % % coef computation
    eq_coefs = zeros(M,3);
    ro = .5;
    for i=1:M
        for ii=1:num_frames
            EQi = 1/ch_resp_est(i,ii);
            EQ_min = 1/ch_resp_est(mod(i-1-1,M)+1,ii);
            EQ_plu = 1/ch_resp_est(mod(i,M)+1,ii);
            eqs = [EQ_min EQi EQ_plu];
   
            if eq_select == 2
                %   % geometric interpolation proposed by Aurelio & Bellanger
                %   % the coefficient computation from same paper
                %   % it's also approach 2 section 4.1.2 from D3.1
                EQ1 = EQ_min*(EQi/EQ_min)^ro;
                EQ2 = EQi*(EQ_plu/EQi)^ro;

            elseif eq_select == 3
                %   % approach 1 D3.1 section 4.1.1 linear interpolation
                EQ1 = (EQ_min+EQi)/2;
                EQ2 = (EQi+EQ_plu)/2;            
            end
            
            eq_coefs(i,1)= ((-1)^(i-1))*((EQ1-2*EQi+EQ2)+j*(EQ2-EQ1))/4;
            eq_coefs(i,2)= (EQ1+EQ2)/2;
            eq_coefs(i,3)= ((-1)^(i-1))*((EQ1-2*EQi+EQ2)-j*(EQ2-EQ1))/4; 
            
            % %     re im re im ...
            equalizer_output =conv(eq_coefs(i,:),scp_data(i,:));
            sp_output(i,1+(ii-1)*2*syms_per_frame:2*syms_per_frame+(ii-1)*2*syms_per_frame)...
                = equalizer_output(2*K+1+(ii-1)*2*(syms_per_frame+1)+1:1+2*K+2*syms_per_frame+(ii-1)*2*(syms_per_frame+1));
            % 2*K+1+1+(ii-1)*2*(syms_per_frame+1):2*K+2*syms_per_frame+1+(ii-1)*2*(syms_per_frame+1)
            % sp_output(i,:) = equalizer_output(i,2*K+1+1:2*K+2*num_symbols+1);
            
        end
    end
elseif eq_select == 4 %no equalizer
    for i=2*K+1:2*K+2*num_symbols
        sp_output(:,i-2*K) = scp_data(:,i);
    end
end