function [ in_range, mx_lvl] = bfs_connectivity(node_tag, Sink_neighbor, each_side)
% ttt = Sink_neighbor
global Node
global k
cnt = 0;
ind = 1;
mx_lvl = 1;
nodes = intersect(Sink_neighbor, node_tag);
if isempty(nodes)
    in_range = 0;
    return 
end
vis = [];
level = [];
for i = 1:200
    vis = [vis 0];
    level = [level 0];
    Node(i).npar = 0;
end
for i = 1:size(nodes, 2)
    vis(nodes(i)) = 1;
    level(nodes(i)) = 1;
    Node(nodes(i)).lvl = 1;
end
% coutNeigh = nodes
while ind <= size(nodes, 2)
    cur = nodes(ind);
    coutcur = cur;
    tt = Node(cur).neighbor
    for i = 1:size(Node(cur).neighbor, 2)
        go = Node(cur).neighbor(i);
        coutgo = go;
        coutvisgo = vis(go);
%         pause;
        if vis(go) == 0 && Node(go).status == 1 && Node(go).exist==1
            vis(go) = 1;
            nodes = [nodes go];
            Node(go).lvl = Node(cur).lvl + 1;
            mx_lvl = max(mx_lvl, Node(go).lvl);
            Node(go).npar = 1;
        else if vis(go) == 1 && Node(go).lvl == Node(cur).lvl + 1
                Node(go).npar = Node(go).npar + 1;
        end
        end
    end
    ind = ind + 1;
end
for i = 1:size(nodes, 2)
    level(nodes(i)) = Node(nodes(i)).lvl;
end
% tt = nodes
% pause;
% for i = 1:k
%     fprintf('%i %i \n',i, Node(i).lvl);
% end
% pause;
for i = size(nodes, 2):-1:1
    cur = nodes(i);
    Node(cur).lvl = 1;
    for cc = 1:size(Node(cur).neighbor, 2)
        go = Node(cur).neighbor(cc);
        if vis(go) == 1 && level(go) == level(cur) + 1
            Node(cur).lvl = Node(cur).lvl + (Node(go).lvl / Node(go).npar);
        end
    end
end
% for i = 1:size(nodes, 2)
%     Node(nodes(i)).lvl = mx_lvl - Node(nodes(i)).lvl + 1;
%     fprintf('%i %i\n', i, Node(nodes(i)).npar);
% end
% pause;
in_range = size(nodes, 2);
end