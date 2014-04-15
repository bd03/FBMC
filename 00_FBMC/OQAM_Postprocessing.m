%% OQAM_Postprocessing
%
% Burak Dayi
%
% This script will perform OQAM postprocessing as it was described in
% PHYDYAS Deliverable 5.1
%
% Dependencies: sp_output - processed signal output 
% Output: oqam_demod - estimated message
%
% Created: 02-03-2014

% disp('OQAM postprocessing')

oqam_demod = zeros(M,2*num_symbols); %the matrix that will store demodulated signal
oqam_input = zeros(M,2*num_symbols); %this will hold after taking real part

mjn= (-j).^((1:2*num_symbols)-1); % (-j)^n

for k=1:M
    
    theta = mjn*(-j)^(k-1);
    oqam_input(k,:) =  real(sp_output(k,:).*theta);
    if ~mod((k-1),2) % even k
%         oqam_input(k,:) = real(sp_output(k,:).*theta);
        oqam_demod(k,:) = oqam_input(k,:).*repmat([1 j],1,num_symbols);
    else % odd k
%         oqam_input(k,:) = real(sp_output(k,:).*theta); 
        oqam_demod(k,:) = oqam_input(k,:).*repmat([j 1],1,num_symbols);
    end
end