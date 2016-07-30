function [sp_output] = func_Subchannel_Processing(rx_output,M,K,num_frames,syms_per_frame,lower_boundary,upper_boundary,sp_params,num_users)
%% func_Analysis_Filter_Bank
%
% Burak Dayi
%
% This function will implement FBMC RX Subchannel Processing 
%
% Created: 24-11-2014

% first check if subchannel processing parameters are all given.
try
    eq_select = sp_params.eq_select;
    preamble_enabled = sp_params.preamble_enabled;
    preamble = sp_params.preamble;
    preamble_sel = sp_params.preamble_sel;
    extra_zero = sp_params.extra_zero;
    zero_pads = sp_params.zero_pads;
    length_preamble = sp_params.length_preamble;
    est_col = sp_params.estimation_column;
catch
    % missing element
    error('There are elements missing in sp_params parameter structure.');
end

%% Subchannel_processing
sumfactor = 0.6863;

if (eq_select == 2) | (eq_select == 3) % three tap equalizer needs help from neighbor subchannels.
    low = lower_boundary-1; % lower boundary
    upp = upper_boundary+1; % upper boundary
else
    low = lower_boundary; % lower boundary
    upp = upper_boundary; % upper boundary
end

scp_input = rx_output(low:upp,2*K-1:end-(2*K-1)+1); % (2K-1)-1 delay removed

if preamble_enabled

    syms_preamble = length_preamble/M;
    allocated_subchannels = upper_boundary-lower_boundary+1;

    % preamble selection --> estimation
    % equalizer selection --> estimation & equalization

    % find estimation column:        

    % separate preamble and data
    for i=1:num_frames
        % i
        scp_preamble(:,i) = scp_input(:,est_col+(i-1)*syms_per_frame*2+(i-1)*syms_preamble); % collects central preamble information
        scp_data(:,1+(i-1)*(syms_per_frame)*2:i*(syms_per_frame)*2) = scp_input(:,1+i*syms_preamble+(i-1)*(syms_per_frame)*2:i*syms_preamble+i*(syms_per_frame)*2); % collects data info
        
        % get estimation
        % preamble(low:upp,est_col)
        ch_resp_est(:,i) = scp_preamble(:,i)./(preamble(low:upp,est_col)*num_users*sumfactor);%*sumfactor);        
        % scp_preamble(:,i)./(preamble(low:upp,est_col))
    end
    % mean(abs(ch_resp_est.*ch_resp_est))
    % save(strcat('ch_resp_est',int2str(num_users),'.mat'),'ch_resp_est')
    % save(strcat('ch_resp_est',int2str(modulation),'.mat'),'ch_resp_est')
    % save(strcat('scp_preamble',int2str(modulation),'.mat'),'scp_preamble')
    % save(strcat('scp_data',int2str(modulation),'.mat'),'scp_data')

    % ch_resp_est
    
    eq_coefs = NaN;
    % get equalizer coefficients & perform equalization
    if eq_select == 1 % 1: one tap
        % one tap equalizer
        eq_coefs=[];
        for ii=1:num_frames
            eq_coefs =[eq_coefs repmat(ch_resp_est(:,ii),1,syms_per_frame*2)];
            save('eq_coefs.mat','eq_coefs')
        end

        % equalization:
        sp_output = scp_data./eq_coefs;
    elseif (eq_select == 2) | (eq_select == 3) %2,3: % three taps
        eq_coefs = zeros(allocated_subchannels,3);
        ro = .5;
        for i=2:allocated_subchannels+1 %1st and last rows are neighbors
            for ii=1:num_frames
                EQi = 1/ch_resp_est(i,ii);
                EQ_min = 1/ch_resp_est(i-1,ii);
                EQ_plu = 1/ch_resp_est(i+1,ii);
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

                eq_coefs(i-1,1)= ((-1)^(i-1))*((EQ1-2*EQi+EQ2)+j*(EQ2-EQ1))/4;
                eq_coefs(i-1,2)= (EQ1+EQ2)/2;
                eq_coefs(i-1,3)= ((-1)^(i-1))*((EQ1-2*EQi+EQ2)-j*(EQ2-EQ1))/4; 

                % %     re im re im ...
                % equalization:
                % equalizer_output =conv(eq_coefs(i,:),rx_output(lower_boundary+(i-1),:));
                % i
                % ii
                % (1+(ii-1)*2*syms_per_frame:ii*2*syms_per_frame)
                % (1+2+(ii-1)*2*(syms_per_frame)+ii*syms_preamble:2+ii*2*(syms_per_frame)+ii*syms_preamble)
                equalizer_output =conv(eq_coefs(i-1,:),scp_input(i,:));
                % size(equalizer_output)
                % size(eq_coefs)
                % size(scp_input)
                sp_output(i-1,1+(ii-1)*2*syms_per_frame:ii*2*syms_per_frame)...
                    = equalizer_output(2+(ii-1)*2*(syms_per_frame)+ii*syms_preamble:1+ii*2*(syms_per_frame)+ii*syms_preamble);
                    % = equalizer_output(1+2+(ii-1)*2*(syms_per_frame)+ii*syms_preamble:2+ii*2*(syms_per_frame)+ii*syms_preamble);
                    % = equalizer_output(2*K+2+1+(ii-1)*2*(syms_per_frame+1)+(ii-1):2*K+2+1+(ii-1)*2*(syms_per_frame+1)+(ii-1)+2*syms_per_frame-1);
            end
        end
    elseif eq_select == 4 % 4: no equalizer
        sp_output = scp_data;
    else
        error('Invalid eq_select')
    end        
else
    % no need to process anything
    % warning(['As no preamble intoduced, equalization was not applied.']);
    sp_output = scp_input;
    % size(sp_output)
end

% disp('+Subchannel processing is done.');