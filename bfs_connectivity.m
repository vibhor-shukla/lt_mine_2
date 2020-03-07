function [ in_range, mx_lvl] = bfs_connectivity(node_tag, Sink_neighbor, each_side)
ttt = Sink_neighbor
global Node
cnt = 0;
ind = 1;
mx_lvl = 1;
nodes = intersect(Sink_neighbor, node_tag);
if isempty(nodes)
    in_range = 0;
    return 
end
vis = [];
for i = 1:each_side * 3;
    vis = [vis 0];
end
for i = 1:size(nodes, 2)
    vis(nodes(i)) = 1;
    Node(nodes(i)).lvl = 1;
end
while ind <= size(nodes, 2)
    cur = nodes(ind);
    for i = size(Node(cur).neighbor, 2)
        go = Node(cur).neighbor(i);
        if vis(go) == 1
            continue;
        else if Node(:).status== 1 && Node(:).exist==1
            vis(go) = 1;
            nodes = [nodes go];
            Node(go).lvl = Node(cur).lvl + 1;
            mx_lvl = max(mx_lvl, Node(go).lvl);
            end
        end
    end
    ind = ind + 1;
end
for i = 1:size(nodes, 2)
    Node(nodes(i)).lvl = mx_lvl - Node(nodes(i)).lvl + 1;
end
in_range = size(nodes, 2);
end