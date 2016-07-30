function out_matrix = subcarrier_demapping(in_matrix,Q,N,start,mapping)
    
    out_matrix = zeros(N,size(in_matrix,2));
    
    if mapping==0 % LFDMA
        out_matrix(:,:)=in_matrix(start:(start+N-1),:);
    elseif mapping==1 % DFDMA
%         out_matrix(offset:Q:size(in_matrix,1),:) = in_matrix(:,:);
        for i=1:N
            out_matrix(i,:) = in_matrix(mod(start-1+(i-1)*Q,size(in_matrix,1))+1,:);
        end
    else
        error('Unsupported mapping');
    end
end