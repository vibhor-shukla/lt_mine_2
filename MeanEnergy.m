function [mean_energy] = MeanEnergy(k)
global Node
num = 0;
den = 0;
for i=1:k-1
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
fid = fopen('test.txt', 'a+');
fprintf(fid, '%d %d\n', num, den);
fclose(fid);
% pause;
end