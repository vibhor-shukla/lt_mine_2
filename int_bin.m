function [ bin ] = int_bin( x,n )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
bin = zeros([1 n]);
i=1;
while x>0
   bin(i) = mod(x,2);
    x=floor(x/2);
   i = i+1;
end
% while i < n
%     bin(i) = 0;
%     i = i+1;
% end
bin = flip(bin,2);
end
