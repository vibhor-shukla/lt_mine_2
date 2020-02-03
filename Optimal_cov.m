function [ opt_cov,index,value ] = Optimal_cov( Total_cover,Load,mode )
%  Finding the optimal cover from given set of connected covers
% index is passed for debugging when needed
value = []; % storing every parameter value for cover set.
opt_cov =[];
index=0;
n = size(Total_cover,2);
global Node
tl = nnz([Node(:).status]&[Node(:).exist]); %total live nodes
min_v = -inf; %check for current maximum value
% min_s = tl; for size constraint

for i=1:n
    temp_c = Total_cover{i}; %temp cover
    % old load value - comment when not used
    temp_l = Load{i}; % temp load
    %     % proceed only for cover set where all nodes are alive and exists.
    %     if nnz([Node(temp_c).status]&[Node(temp_c).exist])==length(temp_c)
    %         %% reassignning value to load based on current live status of nodes.
    %         temp = {Node([Node(temp_c).tag]).neighbor};
    %             % neighbor of each cover vertex
    %             m = size(temp,2);
    %             sz = zeros([1 m]); %storing load on each cover node for current set.
    %             for j=1:m
    %                 sz(j) = nnz([Node(temp{j}).status]&[Node(temp{j}).exist]);
    %             end
    %         temp_l = sz;
    % commented because Optimal_cov is called in each iteration
    %% Parameter calculation
    spare = temp_c(logical(temp_l==0)); % storing load =0 elements
    temp_c(logical(temp_l==0))=[]; % removing load=0 element
    temp_l(logical(temp_l==0))=[];% removing load=0 element
    if(numel(temp_c)~=0)  % checking division by zero
        %current parameter value-> val
        val = sum([Node(temp_c).energy]./(log2(temp_l)+1));  % don't penalize for size one cover
        val = val*log2(tl/numel(temp_c));
        value = [value val];
        %% unused formulas
        % val = sum([Node(temp_c).energy])./sum(temp_l); % may need to change formula
        % val = sum([Node(temp_c).energy]);
        %val = val/(size(temp_c,2)*sum(temp_l));
        %     val = val*tl*sum(temp_l)/size(temp_c,2);
        %  val = val/log(size(temp_c,2));
    else
        val = 0;
        value = [value val];
    end
    if val>=min_v %&&size(temp_c,2)<=min_s
        min_v = val; %updating minval
        % min_s = size(temp_c,2);
        % checking condition when a node have zero load
        if ~isempty(spare)&&~isempty(temp_c)&&connected(Node(temp_c),mode)==1
            %updating optimal cover and corresponding index
            opt_cov = temp_c; %temp_c in minimal set
            index =i;
        else
            %updating optimal cover and corresponding index
            opt_cov = [temp_c spare]; %both temp_c and spare are required for connection
            index =i;
        end
    end
    %     else
    %         value = [value 0]; % for cover set not concidered
    %     end
end

end
