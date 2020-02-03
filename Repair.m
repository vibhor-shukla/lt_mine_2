function [dead] = Repair( energy )

global Node
% repair node which are dead
dead = find(((~[Node(:).status])&[Node(:).exist])==1);
temp_r(1:length(dead)) = {1};
[Node(dead).status]= temp_r{:}; %set status =1
temp_r(1:length(dead)) = {'S'};
[Node(dead).role]= temp_r{:}; % set role =Sensor
temp_r(1:length(dead)) = {energy};
[Node(dead).energy]= temp_r{:}; % set energy=max;

end
