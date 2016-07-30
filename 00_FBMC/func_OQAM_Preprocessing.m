function [oqam_m,num_oqam_subsymbols] = func_OQAM_Preprocessing(qam_m,M,num_frames,low_boundary,high_boundary,preamble_enabled,preamble,syms_per_frame)
%% func_OQAM_Preprocessing
%
% Burak Dayi
%
% This function will implement OQAM preprocessing
%
% Created: 24-11-2014

%% OQAM Preprocessing
% disp('OQAM preprocessing')

num_symbols = num_frames*syms_per_frame; % total number of data symbols

jn = (j.^((1:2*num_symbols)-1)).'; % j to the power of n
%for odd subchannels flipped version of this will be used.

oqam_m = zeros(M, 2*num_symbols); %the matrix that will store modulated signal

symbol_matrix_index = 1; % keeps the index number in qam_m matrix
for k=low_boundary:high_boundary
    theta = jn*(j^(k-1)); %theta multiplier
    if ~mod((k-1),2) % even k --> index spans from 1 to 256
        real_parts = upsample(real(qam_m(symbol_matrix_index,:).'),2);
        imag_parts = circshift(upsample(imag(qam_m(symbol_matrix_index,:).'),2),1);        
    else % odd k
        real_parts = circshift(upsample(real(qam_m(symbol_matrix_index,:)).',2),1);
        imag_parts = upsample(imag(qam_m(symbol_matrix_index,:).'),2);           
    end
    oqam_m(k,:) = (real_parts+imag_parts).*theta;
    symbol_matrix_index = symbol_matrix_index +1;
end

% jp = zeros(size(oqam_m));

if preamble_enabled
    % preamble insertion
    for i=num_frames:-1:1
        oqam_m = [oqam_m(:,1:2*(i-1)*syms_per_frame) preamble oqam_m(:,2*(i-1)*syms_per_frame+1:end)];
%         jp = [jp(:,1:2*(i-1)*syms_per_frame) preamble jp(:,2*(i-1)*syms_per_frame+1:end)];
    end
end
% 
% %%!!!!!!!!!!!!!!!!hack!!!!!!!!!!!!!!!!!!!!!!!!!!
%num_symbols = num_symbols+num_frames*2; %num_frames= no of preambles
% num_oqam_subsymbols = 2*num_symbols+3*num_frames; %num_frames= no of preambles
num_oqam_subsymbols = size(oqam_m,2);
%in this file all the computation will depend on num_oqam_subsymbols.

% %soru: is preamble oqam modulated or to be modulated???????????????
% %ans: the preamble takes 3 subsymbol duration with given setup. that
% %indicates that the preamble will be appended to the oqam modulated data