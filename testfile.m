% if strcmp(Node(5).pos,'roof')
%     disp('0') 
% else
%     disp('1') 
% end

for i=1:100
    Node(i).mode = 1;
end
% x = [1,2];
% y = [3,4];
% z = [5,6];
% d = zeros([1 2^15]);
% for i=1:2^15
%     d(i) = dist(x,y,z);
% end

% check = zeros([1 n]);
% for i = 1:n
% check(i) = size(Prob{i},2);
% end

% for K = t_ins
%     disp(K)
% end
% ch = [1,0,0,0,1;0,0,1,1,0];
% nnz(ch);
% 
% for i = 1:5
%     i
% end
%ch = zeros();
% for i = 1:10
%    [ch(i,:),ch1(i,:)] = datasample(Prob{2000},S_sample,'Weights',Prob{2000});
% end
% temp = load('data3.mat',B2);

% temp = [1,4,6,2,3];
% for k = temp
%     k
% end
% for j = 2:-1:1
%     for k = 1:j*10
%      js(k) = j;
%     end
% end
% 
% T = find(ch(1,:)==1);
% 
% if 0 || 0
%     a = 1
% else
%     a = 0
%     
% end
% prompt = 'Enter the value';
% x = input(prompt)
% i = 1;
% A = 0;
% count = 0;
% while i < 12
%     if mod(i,2) == 0
%         i = i+1;
%         continue;
%     end
%     A = A + i;
%     i = i+1;
%     count = count + 1;
% end
% count

% for i = 1:5
%     for j = i+1:5
%         i + j
%     end
% end


% A = zeros([400 400]);
% A = uint8(A);
% for i = 1:4001
%     X(:,:,i) = A;
% end


% X = 1:10;
% %s = RandStream('mlfg6331_64');
% for i = 1:10
%     %s.Substream = i;
%     [Z2(i,:),dk2(i,:)] = datasample(X,5,'Replace',true);
% end
% %[Z,dk] = datasample(s,X,5);
%%

%P = zeros([1 2]);
% for i = 1:exp(8)
%     for j = 1:exp(8)
%         P(i,j) = i+j;
%     end
% end
% function [A,n] = testfile(x,y)
% A{1} = x+y;
% A{2} = x.*y;
% n = 4;
% end
%%

%[X,Y] = meshgrid(0:.5:4,-2:.5:2);
% [X_1,Y_1] = meshgrid(X,Y);
        %domain_2 = 0:.03:15;
        %IC1 = [X,Y];
        % fp = fopen('Adjacency.txt','w');
        % fprintf(fp,'%d\n',A);
        % fclose(fp);
        
%% from stream
%a1 = b .* sin(k.*x);
%d1 = (k^2) .* (b.^2) .* cos(k.*x) .* cos(k.*x);
%Y = -tanh((y-a1)./((1+d1).^(1/2))) + c.*y;

%% from lcc
%     display(x);
%can be node+1
        %% testing inc
% xticks(0:15:P(2,size(P,2)));
% yticks(0:50:P(1,size(P,2)));
% set(gca,'YTick',0:50:P(1,size(P,2)));
%% from dist_adj
        % fp = fopen('test2.txt','w');
% fprintf(fp,'%d\n',[D]);
% fclose(fp)
        % fp = fopen('test1.txt','w');
% fprintf(fp,'%d\n',D);
% fclose(fp);
%% plot locations
% figure;
% scatter(x_1,y_1,1,[0 1 0]);
% axis([0 80 -4 4]);
%% testing plot
% stem(x_1,y_1,'b');
% figure;
% scatter(X,Y,1);
% axis([0 100 -4 4]);
% figure;
% I = 1;
% while (I < size(t)+1)
%     T = t(I);
%    % X = x(I);
%    % Y = y(I);
%     ez(@(x,y)stream(x,y,T),[x,y]);
%     I = I+1;
%end

%% stream plot
% X = 0:0.1:70
% Y = stream

%% Testing differentiation
%syms x y t k c B
% f = -tanh((y - B*sin(k*x - c*k*t))/(B^2*k^2*cos(k*x - c*k*t)^2 + 1)^(1/2));
% k = 2*pi/7.5;
% c = 0.12;
% B = B_mod(0);
% t = 10;
%f =@(x,y) -tanh((y - B.*sin(k.*x))./((B.^2).*(k.^2)*(cos(k.*x)).^2 + 1).^(1/2)) + c.*y;
% fimplicit(f,[0 70 -4 4],'b')
% g = diff(f,x)
%g = -(tanh((y - B*sin(k*x - c*k*t))/(B^2*k^2*cos(k*x - c*k*t)^2 + 1)^(1/2))^2 - 1)*((B*k*cos(k*x - c*k*t))/(B^2*k^2*cos(k*x - c*k*t)^2 + 1)^(1/2) - (B^2*k^3*cos(k*x - c*k*t)*sin(k*x - c*k*t)*(y - B*sin(k*x - c*k*t)))/(B^2*k^2*cos(k*x - c*k*t)^2 + 1)^(3/2));
% h = diff(f,y);
%h = (tanh((y - B*sin(k*x - c*k*t))/(B^2*k^2*cos(k*x - c*k*t)^2 + 1)^(1/2))^2 - 1)/(B^2*k^2*cos(k*x - c*k*t)^2 + 1)^(1/2);


%syms x y t k c b
%t = 1.5;
%c = 0.12;
% f = tanh((y-c)./t);
% g = diff(f,y);
% g
% vpa(subs(g,y,2))
% y = 2;
% h = -(tanh((c - y)/t)^2 - 1)/t
% check the square in h;
% f = tanh((y - b*sin(k*(x - c*t)))/((1 + b*cos(k*(x - c*t))*cos(k*(x - c*t)))^(1/2)));
% g = diff(f,x);
% g
% h = (tanh((y - b*sin(k*(x - c*t)))/(b*cos(k*(x - c*t))^2 + 1)^(1/2))^2 - 1)*((b*k*cos(k*(x - c*t)))/(b*cos(k*(x - c*t))^2 + 1)^(1/2) - (b*k*cos(k*(x - c*t))*sin(k*(x - c*t))*(y - b*sin(k*(x - c*t))))/(b*cos(k*(x - c*t))^2 + 1)^(3/2));

%x = [2,23;43,3;90,3;564,80;1,2000;9,4212]
%max(x(3,:))

%% extra from Are_rec
%N = [];
%N = [N temp1 0];
% N = [N -1];
% x_range = x_min:x_count:x_max;
% y_range = y_min:y_count:y_max;
% x_range = x_min.:(x_max-x_min)./div.:x_max;
% y_range = y_min.:(y_max-y_min)./div.:y_max;
% n = size(x_range,2);
% m = size(y_range,2);
%   P = zeros([1 n-1*m-1]); %cell();
%% from Dist_adj2
% DO IT OUTSIDE
    %D(:,:,k) = D(:,:,k) + transpose(D(:,:,k)) - diag(diag(D(:,:,k)));
    %% extra edge_red
    % do this outside
   % R_1 = R_1 + transpose(R_1) - diag(diag(R_1));
   %     R_1 = ones([b b],'uint8');
   %global pos
    
%% from edgecount_1

%     if count{I}(k) == 0
%         P(k,:,2) = k-1;
%         P(k,k,2) = 0;
%     end

%% from edge_count

% W = zeros([1 n]);
% X = zeros([1 n]);
% Y = zeros([1 n]);
% Z = zeros([1 n]);
%global pos
%% extra from Algo
% size exceding have to store them seperately. and access individually.
% A = Dist_adj1(x_1,y_1);
% A(:,:,k)= A(:,:,k) - diag(diag(A(:,:,k)));

% work need to be done.
% currently using uint8 to reduce space in Dist_adj1 function
%[Prob,Loc] = Area_rec(x_1,y_1);
% Global has less time
% B = Edge_1(A,sz);  first case%%  pass pos in it to reduce time, but it increases
% A(:,:,k)= A(:,:,k) - diag(diag(A(:,:,k)));
%s_2(K) = RandStream('mlfg6331_64'); %use substream in loop
%t_index = K;
%t_ins = domain_1(K);
% s_2(K).Substream = I;
%Lcc_store(c) = Lcc;
%Lcc_store(c) = Lcc;
%alpha = avg_Lcc1-avg_Lcc2;
%hold;
%plot();
%P = 1:200; %experimenting
%plot(avg_Lcc1,three_1,'r');
% figure;
% P = t_ins;
% plot(domain_1(P),avg_Lcc1(P),'r',domain_1(P),avg_Lcc2(P),'b');
%figure();
%plot(domain_1(P),alpha,'g');
%plot(domain_1(P),new1(P),'r',domain_1(P),new2(P),'b');
%plot(domain_1(P),lcc20(P),'g',domain_1(P),lcc60(P),'r',domain_1(P),lcc80(P),'b',domain_1(P),avg_Lcc2(P),'m');
%plot(domain_1(P),avg_Lcc1(P),'r',domain_1(P),hh(P),'b');
%plot(domain_1(P),avg_Lcc1(P),'r',domain_1(P),up(P),'b',domain_1(P),hh(P),'g',domain_1(P),low(P),'m');

%% Algorithmic process for single database entry and one run of probabilistic model.

% Create a uniform/normal-random/defined inital distrbution using random number generator.

% Run the mendering model to get the final distribution at period of interest and store in a file/database/variable

% Apply the area division algorithm and generate the locality set and calculate the probability for them.

% Choose the network state from the possibilities by efficient sampling or sir method of choosing a random number with a constraints.

% Use the edge relation to calculate the probability of connectivity and compare it with real connectivity.
% Comparing has to incorporate all the possible location utilization in calculation.

%% extra from metric

% for I = 1:size(B1,2)
%    % temp = cast(B1{I},'single');
%    % e1{I} = eig(temp);
%    e1{I} = eig(cast(B1{I},'single'));
% end
% 
% for I = 1:size(B2,2)
% %     temp = cast(B2{I},'single');
% %     e2{I} = eig(temp);
%     e2{I} = eig(cast(B2{I},'single'));
% end
% 
% for I = 1:size(B3,2)
% %     temp = cast(B3{I},'single');
% %     e3{I} = eig(temp);
%     e3{I} = eig(cast(B3{I},'single'));
% end

% clear temp