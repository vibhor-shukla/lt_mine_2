% function [  ] = Energy_fun( node_tag,t,it )
%Update Energy of live nodes having tag(node_tag) at time(t) and after (it) iterations.
% function [  ] = Energy_fun( node_tag )
function [] = Energy_fun(flag, each_count, mx_lvl)
%flag to indicate the instance of called state
global Node
node_tag = find([Node(:).exist]&[Node(:).status]==1);
% Time variation is not considered as energy consumption is sumilated as a
% fixed amount for each work type.
S_fxd = 1.047745358e-3;
R_fxd = 2.166226913e-3;
E_fxd = 3.5e-3;
D_fxd = 4e-3;
E_threshold = 2.5;
E_die = 1.8;
n = length(node_tag);
reduce = (S_fxd + R_fxd);
if flag==1 %  when health state called
    for i=1:n
        if (Node(node_tag(i)).energy<=E_die) % condition for dead node
            Node(node_tag(i)).status=0;
            continue;
        end
        Node(node_tag(i)).energy = Node(node_tag(i)).energy-D_fxd ;
        if (Node(node_tag(i)).energy<=E_die) % condition for dead node
            Node(node_tag(i)).status=0;
        end
    end
% else % Normal state
%     for i=1:n %for all live node
%         if (Node(node_tag(i)).energy<=E_die) % condition for dead node
%             Node(node_tag(i)).status=0;
%             continue;
%         end
%         % checking for Role conversion
% %         if ((Node(node_tag(i)).energy<=E_threshold) && (Node(node_tag(i)).role~='S'))
% %             Node(node_tag(i)).role='S';
% %         end
%         if Node(node_tag(i)).role=='S'% for sensors
%             Node(node_tag(i)).energy = Node(node_tag(i)).energy-S_fxd ;
%             if (Node(node_tag(i)).energy<=E_die) % condition for dead node
%                 Node(node_tag(i)).status=0;
%             end
%         elseif Node(node_tag(i)).role=='R'%for relays
%             Node(node_tag(i)).energy = Node(node_tag(i)).energy-R_fxd ;
%             if (Node(node_tag(i)).energy<=E_die) % condition for dead node
%                 Node(node_tag(i)).status=0;
%             end
%         elseif Node(node_tag(i)).role=='E'%for emergency
%             Node(node_tag(i)).energy = Node(node_tag(i)).energy-E_fxd;
%             if (Node(node_tag(i)).energy<=E_die) % condition for dead node
%                 Node(node_tag(i)).status=0;
%             end
%         elseif Node(node_tag(i)).role=='B'%for both
%             Node(node_tag(i)).energy = Node(node_tag(i)).energy-S_fxd-R_fxd;
%             if (Node(node_tag(i)).energy<=E_die) % condition for dead node
%                 Node(node_tag(i)).status=0;
%             end
%         end
%     end
% end
else
    for i = 1:n
        if (Node(node_tag(i)).energy<=E_die) % condition for dead node
             Node(node_tag(i)).status=0;
        else
%              Node(node_tag(i)).energy = Node(node_tag(i)).energy-S_fxd-R_fxd;
            Node(node_tag(i)).energy = Node(node_tag(i)).energy-(reduce * Node(node_tag(i)).lvl);
            if (Node(node_tag(i)).energy<=E_die) % condition for dead node
                Node(node_tag(i)).status=0;
            end
        end
    end
end
