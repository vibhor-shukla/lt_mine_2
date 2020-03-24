function [ vertex,load ] = cover( l, node_tag,Sink_neighbor,mode )
% Find all conncted vertex cover

global Node
vertex = {}; % output cell of possible cover sets
load = {}; % load of individual node in cover set
n = size(node_tag,2); % total live nodes
sum1 = sum(node_tag); % Sum of all tag values is unique.

if n == 0
    return;
end

% For all possible power set of live nodes being as cover.
%for i = 1:(2^n-1);
%  i = (2^n - 1); 
    % decimal to binary
%     pos_cov = de2bi(i,n); %LSB-MSB
    pos_cov = [];
    for i = 1:n
        pos_cov = [pos_cov 1];
    end
    %pos_cov
    %pos_cov = int_bin(i,n); MSB-LSB
    pos_cov = logical(pos_cov);
    % collection of all neighbors of taken cover set.
    pos_neibr = [Node(node_tag(pos_cov)).neighbor];
    temp = unique(pos_neibr); % finding unique neighbors for all possible set
    temp = temp(logical([Node(temp).status]&[Node(temp).exist]==1)); % only live neighbor
    if size(temp,2)==size(node_tag,2) && sum(temp)==sum1 % only when covering all nodes sum is equal otherwise less
        %only cover which have direct commuication to sink
        if numel(intersect(Sink_neighbor,[Node(node_tag(pos_cov)).tag]))>0
            % conditioned only for connected cover
            %conn = Node(node_tag(pos_cov)); % cover vertices
            conn = connected(Node(node_tag(pos_cov)),mode);
            if conn==1
                vertex{end+1} = [Node(node_tag(pos_cov)).tag]; %cover set
                temp = {Node([Node(node_tag(pos_cov)).tag]).neighbor};
                % neighbor of each cover vertex
                m = size(temp,2);
                sz = zeros([1 m]); %storing load on each cover node for current set.
                for j=1:m
                    sz(j) = nnz([Node(temp{j}).status]&[Node(temp{j}).exist]);
                end
                load{end+1} = sz; % load of nodes in cover set.
            end
        end
    end
%end
end
