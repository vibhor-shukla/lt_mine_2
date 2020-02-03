function [Total_node, X, Y, Z,each_side,Sink] = Deploy(length_interest,Rcom)
%% input variables
% clear
%length_interest = 100;   
%length of deployment corridor. 
%global Rcom
%Rcom=40; % communication  range.
Height=3.8; % height of tunnel.
h3 = Height / 3;
Width=5;   % width of tunnel.
%x_int = 90*Rcom/100; % deployment distance of nodes in length.
%x_int = sqrt(((0.8 * Rcom) ^ 2) - (Width * Width));
x_int = 0.45 * Rcom;
%% Deciding coordinates of nodes.
% w1-> wall 1,  w2->wall 2, r->roof.

X_w1 = 0:2 * x_int:length_interest;
%X_r = X_w1(1)+(2 * x_int / 3):2 * x_int:length_interest-x_int/2;
X_w2 = x_int:2 * x_int:length_interest;

s1 = size(X_w1, 2);
%s2 = size(X_r, 2);
s3 = size(X_w2, 2);

%siz = min(min(s1, s2), s3);
siz = min(s1, s3);

while(size(X_w1, 2) > siz)
    X_w1(end) = [];
end

%while(size(X_r, 2) > siz)
%    X_r(end) = [];
%end

while(size(X_w2, 2) > siz)
    X_w2(end) = [];
end

%size(X_w1, 2) + size(X_r, 2) + size(X_w2, 2)
size(X_w1, 2) + size(X_w2, 2)

Y_w1 = zeros([1 size(X_w1,2)]);
%Y_r = (Width/2).*ones([1 size(X_r,2)]);
Y_w2 = Y_w1 + Width;

Z_w1 = (Height / 2).*ones([1 size(X_w1,2)]);
%Z_r = Height.*ones([1 size(X_r,2)]);
Z_w2 = (Height / 2).*ones([1 size(X_w2,2)]);

for i = 1:size(X_w1, 2)
    x = mod(i - 1, 3);
    if (x == 0)
        Z_w1(i) = h3;
    elseif (x == 1)
        Z_w1(i) = h3 * 3;
    else
        Z_w1(i) = h3 * 2;
    end
end

for i = 1:size(X_w2, 2)
    x = mod(i - 1, 3);
    if (x == 0)
        Z_w2(i) = h3 * 2;
    elseif (x == 1)
        Z_w2(i) = h3;
    else
        Z_w2(i) = h3 * 3;
    end
end

% storing all X,Y,Z indices of nodes
X = {X_w1;X_w2};
Y = {Y_w1;Y_w2};
Z = {Z_w1;Z_w2};

% counting number of nodes deployed
n_w1 = size(X_w1,2);
n_w2 = size(X_w2,2);
%n_r = size(X_r,2);
Total_node = n_w1+n_w2;
each_side = [n_w1,n_w2];

% initializing node degrees to zeros.
deg_w1=zeros([1 n_w1]); 
deg_w2=zeros([1 n_w2]);
%deg_r=zeros([1 n_r]);

% sink
X_s = X_w1(end)+(x_int/2);
Y_s = Width / 2;
Z_s = h3 * 2;
Sink = [X_s,Y_s,Z_s];
%% plotting.
%plotting nodes and defining connectivity between them
view(0, 90)
figure(1);
scatter3(X_w1,Y_w1,Z_w1,'b');
hold on;
scatter3(X_w2,Y_w2,Z_w2,'b');
%scatter3(X_r,Y_r,Z_r,'r');
scatter3(X_s,Y_s,Z_s,'black','s','fill');

for i = 1:n_w1
    % betweem wall 1 and wall 2.
    for j= 1:n_w2
        if dist([X_w1(i) X_w2(j)],[Y_w1(i) Y_w2(j)],[Z_w1(i) Z_w2(j)])<=Rcom
             plot3([X_w1(i) X_w2(j)],[Y_w1(i) Y_w2(j)],[Z_w1(i) Z_w2(j)],'b');
             % updating degree count of nodes.
             deg_w1(i)=deg_w1(i)+1;
             deg_w2(j)=deg_w2(j)+1;
        end
    end 
    % Between wall 1 and roof
%     for j= 1:n_r
%         if dist([X_w1(i) X_r(j)],[Y_w1(i) Y_r(j)],[Z_w1(i) Z_r(j)])<=Rcom
%              plot3([X_w1(i) X_r(j)],[Y_w1(i) Y_r(j)],[Z_w1(i) Z_r(j)],'r');
%              % updating degree count of nodes.
%              deg_w1(i)=deg_w1(i)+1;
%              deg_r(j)=deg_r(j)+1;
%         end
%     end
end

% Between wall 2 and roof
% for i= 1:n_w2
%     for j= 1:n_r
%         if dist([X_w2(i) X_r(j)],[Y_w2(i) Y_r(j)],[Z_w2(i) Z_r(j)])<=Rcom
%              plot3([X_w2(i) X_r(j)],[Y_w2(i) Y_r(j)],[Z_w2(i) Z_r(j)],'r');
%              % updating degree count of nodes.
%              deg_w2(i)=deg_w2(i)+1;
%              deg_r(j)=deg_r(j)+1;
%         end
%     end
% end

 xlim([0 length_interest]);
 ylim([0 Width]);
 zlim([0 Height]);

% sphere around sink
[x,y,z]=sphere(128);
h=surfl(x*Rcom+Sink(1),y*Rcom+Sink(2),z*Rcom+Sink(3));
set(h,'FaceAlpha',0.1);
shading interp;
hold off;

% plotting degree of nodes in the graph
figure(2);
subplot(1,3,1)
plot(1:n_w1,deg_w1,'-*');
subplot(1,3,2)
plot(1:n_w2,deg_w2,'-*');
% subplot(1,3,3)
% plot(1:n_r,deg_r,'-*');

%% Rough work testing functions.

%n_x_w = floor(length_interest*2/Rcom + 1);
%n_x_r = floor(length_interest*2/Rcom);
%X_w1=zeros([1 n_x_w]);

% figure;
% 
% X=[1,2,3];
% Y=[1,2,3];
% Z=[1,2,3];
% scatter3(X,Y,Z);
% hold;
% plot3(X,Y,Z);

%axis([0 lenght_interest 0 Width 0 Height]);
% for i = 1:n_x_w
%     X_
% end




end