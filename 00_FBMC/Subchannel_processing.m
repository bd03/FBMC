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

scp_input = rx_output(:,2*K-1:end-(2*K-1)+1);
scp_preamble = zeros(M,num_frames);
scp_data = zeros(M,num_symbols*2);
ch_resp_est =zeros(M,num_frames);
sp_output = zeros(M,num_symbols*2);

if strcmp(estimation_method,'IAM')
    % remove preamble
    for i=1:num_frames
        scp_preamble(:,i) = scp_input(:,2+(i-1)*(syms_per_frame+1)*2+(i-1));
        scp_data(:,1+(i-1)*(syms_per_frame)*2:i*(syms_per_frame)*2) = scp_input(:,4+(i-1)*(syms_per_frame+1)*2+(i-1):4+(i-1)*(syms_per_frame+1)*2+(i-1)+syms_per_frame*2-1);
    end
    
    % estimation of channel
    for i=1:num_frames
        ch_resp_est(:,i) = scp_preamble(:,i)./(preamble(:,2)*sumfactor);
    end

    if eq_select == 1
        % one tap equalizer
        eq_coefs=[];
        for ii=1:num_frames
            eq_coefs =[eq_coefs repmat(ch_resp_est(:,ii),1,syms_per_frame*2)];
        end

        sp_output = scp_data./eq_coefs;
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
                equalizer_output =conv(eq_coefs(i,:),rx_output(i,:));
                sp_output(i,1+(ii-1)*2*syms_per_frame:2*syms_per_frame+(ii-1)*2*syms_per_frame)...
                    = equalizer_output(2*K+2+1+(ii-1)*2*(syms_per_frame+1)+(ii-1):2*K+2+1+(ii-1)*2*(syms_per_frame+1)+(ii-1)+2*syms_per_frame-1);
            end
        end
    elseif eq_select == 4 %no equalizer
        sp_output=scp_data;
    end
elseif strcmp(estimation_method,'IAM4')
    % remove preamble
    for i=1:num_frames
        scp_preamble(:,i) = scp_input(:,3+(i-1)*(syms_per_frame+2)*2);
        scp_data(:,1+(i-1)*(syms_per_frame)*2:i*(syms_per_frame)*2) = scp_input(:,5+(i-1)*(syms_per_frame+2)*2:5+(i-1)*(syms_per_frame+2)*2+syms_per_frame*2-1);
    end
    
    %
    % %estimation of channel
    for i=1:num_frames
        ch_resp_est(:,i) = scp_preamble(:,i)./(preamble(:,3)*sumfactor);
    end

    if eq_select == 1
        % one tap equalizer
        eq_coefs=[];
        for ii=1:num_frames
            eq_coefs =[eq_coefs repmat(ch_resp_est(:,ii),1,syms_per_frame*2)];
        end

        sp_output = scp_data./eq_coefs;
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
                equalizer_output =conv(eq_coefs(i,:),rx_output(i,:));
                sp_output(i,1+(ii-1)*2*syms_per_frame:2*syms_per_frame+(ii-1)*2*syms_per_frame)...
                    = equalizer_output(2*K+2+1+(ii-1)*2*(syms_per_frame+2)+1:1+2+2*K+2*syms_per_frame+(ii-1)*2*(syms_per_frame+2));
            end
        end
    elseif eq_select == 4 %no equalizer
        sp_output=scp_data;
    end
end