function [ ] = Destroy_wall( wall )

global Node
Pos = {'wall_1','roof','wall_2'};
x = Pos(wall); %choosen side
n = length(Node);
for i=1:n
    if ~isempty(Node(i).pos) && Node(i).pos==x && Node(i).exist==1
        Node(i).status = 0;  % allocating all node to dead on choosen side.
    end
end

end
