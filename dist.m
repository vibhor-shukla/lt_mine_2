function [ D ] = dist( X,Y,Z )
% Distance between two coordinates

D=( (X(1)-X(2))^2 + (Y(1)-Y(2))^2 + (Z(1)-Z(2))^2 ).^(1/2);

end
