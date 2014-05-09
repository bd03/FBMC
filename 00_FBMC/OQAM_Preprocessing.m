% OQAM_Preprocessing
%
% Burak Dayi
%
% This script will perform OQAM preprocessing as it was described in
% PHYDYAS Deliverable 5.1
%
% Dependencies: qam_m - QAM modulated message, M, num_symbols
% Output: oqam_m - OQAM modulated message
%
% Created: 02-03-2014

% disp('OQAM preprocessing')

jn = (j.^((1:2*num_symbols)-1)).'; % j to the power of n
%for odd subchannels flipped version of this will be used.

oqam_m = zeros(M, 2*num_symbols); %the matrix that will store modulated signal

for k=1:M
    theta = jn*(j^(k-1)); %theta multiplier
    if ~mod((k-1),2) % even k
        real_parts = upsample(real(qam_m(k,:).'),2);
        imag_parts = circshift(upsample(imag(qam_m(k,:).'),2),1);        
    else % odd k
        real_parts = circshift(upsample(real(qam_m(k,:)).',2),1);
        imag_parts = upsample(imag(qam_m(k,:).'),2);           
    end
    oqam_m(k,:) = (real_parts+imag_parts).*theta; 
end

jp = zeros(size(oqam_m));

% preamble insertion
for i=num_frames:-1:1
    oqam_m = [oqam_m(:,1:2*(i-1)*syms_per_frame) preamble oqam_m(:,2*(i-1)*syms_per_frame+1:end)];
    jp = [jp(:,1:2*(i-1)*syms_per_frame) preamble jp(:,2*(i-1)*syms_per_frame+1:end)];
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