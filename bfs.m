function [ ] = bfs( v,d)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
global visited
visited(v) = 1;
n = size(d,1);
for i=1:n
    if d(v,i)==1&&visited(i)==0
        bfs(i,d);
    end
end
