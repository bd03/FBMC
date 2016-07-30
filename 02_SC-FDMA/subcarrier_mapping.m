function out_matrix = subcarrier_mapping(in_matrix,Q,M,start,mapping)
    
    out_matrix = zeros(M,size(in_matrix,2));
    
    if mapping==0 % LFDMA
        out_matrix(start:(start+size(in_matrix,1)-1),:)=in_matrix(:,:);
    elseif mapping==1 % DFDMA
%         out_matrix(offset:Q:size(in_matrix,1),:) = in_matrix(:,:);
        for i=1:size(in_matrix,1)
            out_matrix(mod(start-1+(i-1)*Q,M)+1,:) = in_matrix(i,:);
        end
    else
        error('Unsupported mapping');
    end
end