function [mean_energy] = MeanEnergy(k)
global Node
num = 0;
den = 0;
for i=1:k
    if Node(i).status == 1 && Node(i).exist == 1
        num = num + Node(i).energy;
        den = den + 1;
    end
end
mean_energy = 0;
if den == 0
    return ;
end
mean_energy = num / k;
end