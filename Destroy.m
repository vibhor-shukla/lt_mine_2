function [ tag ] = Destroy( percentage )
% wall collapse feature not added
global Node
tag =[];
%rng('shuffle');
live = [Node(:).exist]&[Node(:).status];
n = nnz(live); %total live node
% percentage = 40; % maximum allowed percentage of node to destroy.
n = uint16(percentage*n/100);
if(n~=0)
    r = randi(n,1);% random number of node to distroy
    loc = find(live==1); % Live node index .
    %enrgy = [Node(live).energy];
    %randsample(loc,r,false,enrgy); %weighted sampling without replacement is not supported.
    sample = randsample(loc,r,false); %sampled node to destroyed.
    for  i=1:length(sample)
            Node(sample(i)).status = 0;  % allocating all node to dead;
            Node(sample(i)).exist = 1;
    end
    tag = [Node(sample).tag];
end

end
