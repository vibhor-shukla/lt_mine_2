function [  ] = Addition(Length_interest,current,new)
%Input is current length of deployment and additional length 'new' required
global Node
n = length(Node);
% length adjustment.
start = Length_interest-current;
END = start-new;
if start>0 %restricting over calculation when Length_interest==current
    for i=1:n
        x=[Node(i).loc];
        if(numel(x)>0)
            x = x(1);
            %condition for node lying in region Staring from sink side i.e. far end
            if(x<=start&&x>=END)
                Node(i).exist=1;
                Node(i).status=1;
            end
        end
    end
end

end
