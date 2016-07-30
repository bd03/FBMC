function [indices]=determine_indices(num_users,subch_per_user,M,approach)
% this function outputs a vector of indices within [1,M] that specifies the
% boundaries for different users

min_subch = num_users*subch_per_user+num_users;

if M<min_subch
    error('M is smaller than the minimum number of subchannel required');
end

indices = zeros(1,2*num_users);

% two methods can be adopted. one min seperation of the users and have
% maximum separation on the edges. The other approach would be equal
% separation of each user.



if approach == 1
    % approach 1 -ok
    low_nonallocated = floor((M-min_subch)/2);
    high_nonallocated = M-low_nonallocated;

    indices(1,1) = low_nonallocated;

    for i=2:2*num_users
        if mod(i,2)==0
            indices(1,i)=indices(1,i-1)+subch_per_user-1;
        else
            indices(1,i)=indices(1,i-1)+2;
        end
    end

    indices = 1+indices;
elseif approach == 2
    % approach 2
    total_nonallocated = M-num_users*subch_per_user;
    base_separation = floor(total_nonallocated/num_users);
%     first_separation = floor(base_separation/2);
    residual = mod(total_nonallocated,num_users);
    
    indices(1,1) = 1;
    
    if residual>0
        indices(1,1) = indices(1,1)+1;
        residual = residual - 1;
    end
    
    for i=2:2*num_users
        if mod(i,2)==0
            indices(1,i)=indices(1,i-1)+subch_per_user-1;
        else
            indices(1,i)=indices(1,i-1)+base_separation+1;
            if residual>0
                indices(1,i) = indices(1,i)+1;
                residual = residual - 1;
            end
        end
    end
    
    indices = floor(base_separation/2)+indices;
else
    error('approach should either be 1 or 2.');
end






