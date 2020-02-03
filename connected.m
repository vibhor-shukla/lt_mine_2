function [ b ] = connected( conn,mode )
% Check the connectivity of current vertex cover
% input structure of node

n = size(conn,1);
m = size(conn,2);
n = max(n,m);
global Rcom
% storing location of nodes
Loc = zeros([n 3]);
for i=1:n
    Loc(i,:) = conn(i).loc;
end

%Initializing direct connectivity
d = zeros([n n]);
for i=1:n
   for j=1:n
       v = dist([Loc(i,1) Loc(j,1)],[Loc(i,2) Loc(j,2)],[Loc(i,3) Loc(j,3)]);
       if v<=Rcom&&mode %mode=1
            d(i,j)=1;
       elseif v<=Rcom&&~strcmp(conn(i).pos,conn(j).pos) %mode=0
           d(i,j)=1;
       end
   end
end

global visited
visited = zeros([1 n]);
bfs(1,d);
if nnz(visited)==n
    b=1;
else
    b=0;
end
% % transitive closure
% for k=1:n
%     for i=1:n
%         for j=1:n
%             d(i,j) = d(i,j)||(d(i,k)&&d(k,j));
%         end
%     end
% end
% 
% % checking connectivity in cover set
% if nnz(d)==n^2
%     b=1;
% else
%     b=0;
% end

end
