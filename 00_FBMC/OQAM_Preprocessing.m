%% OQAM_Preprocessing
%
% Burak Dayi
%
% This script will perform OQAM preprocessing as it was described in
% PHYDYAS Deliverable 5.1
%
% Dependencies: qam_m - QAM modulated message, M, num_symbols
% Output: oqam_m - OQAM modulated message
%
% Last updated: 18-03-2014

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


