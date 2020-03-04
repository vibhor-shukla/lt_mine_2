function [ in_range] = bfs_connectivity(node_tag, Sink_neighbor, each_side)

global Node
cnt = 0;
ind = 1;
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
    vis[nodes(i)] = 1;
end
while ind <= size(nodes, 2)
    cur = nodes(ind);
    for i = size(Node(cur).neighbor, 2)
        go = Node(cur).neighbor(i);
        if vis(go) == 1
            continue;
        else
            vis(go) = 1;
            nodes = [nodes go];
        end
    end
    ind = ind + 1;
end
in_range = size(nodes, 2);
end