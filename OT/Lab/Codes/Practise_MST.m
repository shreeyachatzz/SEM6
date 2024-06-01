%% GRAPHICAL
%% QUESTION 1: Graphical method to solve 
% Max Z= 6x1+11x2
% 2x1+x2<=104
% x1+2x2<=76
% x2>=0
% x1>=0
clf
clc
clear all
format short

%(1)Describing the question
c=[6,11];
A=[2,1;1,2;0,1;1,0];
B=[104;76;0;0];
n=size(A,1);
const=[1,1];
objective=1;
x1=0:0.01:max(B);

%(2)Finding y values
for i=1:n-2
    y(i,:)=(B(i)-A(i,1)*x1)/A(i,2);
end
%(3)Plotting the lines
for i=1:n-2
    y(i,:)=max(0,y(i,:));
    plot(x1,y(i,:),"LineWidth",4);
    hold on
end
hold on

%(4)Finding the points of intersection
pt=[0;0];
for i=1:n
    A1=A(i,:);
    B1=B(i);
    for j=i+1:n
        A2=A(j,:);
        B2=B(j);
        A3=[A1;A2];
        B3=[B1;B2];
        X3=A3\B3;
        if (X3>=0)
            pt=[pt X3];
        end
    end
end

%(5)Inversing pt and storing it in X and keeping only unique points
X=pt'
X=unique(X,'rows');
hold on

%(6)Finding the feasible points
x1=X(:,1);
x2=X(:,2);
for i=1:n-2
    if const(i)==1
        ind=find(A(i,:)*X'>B(i));
        X(ind, :)=[];
    else
        ind=find(A(i,:)*X'<B(i));
        X(ind, :)=[];
    end
end

%(7) Finding the optimal result
obj_val=c*X';
if(objective==1)
    [value,ind]=max(obj_val);
else
    [value,ind]=min(obj_val);
end
fprintf("Optimal Value: %f\n",value);
fprintf("Optimal Point: %g %g\n",X(ind,:));

%(8) Labelling the graph and plotting it
x=X(:,1);
y=X(:,2);
scatter(x,y,'*');
hold on
k=convhull(x,y);
fill(x(k),y(k),'m');
xlim([0 max(x)+1]);
ylim([0 max(y)+1]);
title("Plot for the LPP");
xlabel('X AXIS ->');
ylabel('Y AXIS ->');

%(9) Checking the solution
x=0:0.1:max(B);
for z=0:8:value
    y=(z-c(1)*x)/c(2);
    plot(x,y);
    hold on
    drawnow
    pause(0.001)
end
hold on

%% QUESTION 2: Graphical method to solve 
% Max Z= 5x1+8x2
% 1x1+2x2<=200
% x1+1x2<=150
% x2<=60
% x1,x2>=0

clf
clc
clear all
format short
%(1) Define the function
c=[5,8];
A=[1,2;1,1;0,1;0,1;1,0];
B=[200;150;60;0;0];
const=[1,1,1];
objective=1;
n=size(A,1);
x1=0:0.01:max(B);

%(2) Find the values of y
for i=1:n-2
    y(i,:)=(B(i)-A(i,1)*x1)/A(i,2);
end

%(3) Plot the lines
for i=1:n-2
    y(i,:)=max(0,y(i,:));
    plot(x1,y,'LineWidth',4);
    hold on
end
hold on

%(4) Find the points of intersection
pt=[0;0];
for i=1:n
    A1=A(i,:);
    B1=B(i);
    for j=i+1:n
        A2=A(j,:);
        B2=B(j);
        A3=[A1;A2];
        B3=[B1;B2];
        X3=A3\B3;
        if X3>=0
            pt=[pt X3];
        end
    end
end

%(5) Inverse pt for X and keep only the unique points
X=pt';
X=unique(X,'rows');

%(6) Keep only the feasible points
for i=1:n-2
    if(const(i)==1)
        ind=find(A(i,:)*X'>B(i));
        X(ind,:)=[];
    else
        ind=find(A(i,:)*X'<B(i));
        X(ind,:)=[];
    end
end

%(7) Find the optimal values
obj_val=c*X';
if objective==1
    [value,ind]=max(obj_val);
else
    [value,ind]=min(obj_val);
end

fprintf("Optimal Value: %f\nOptimal Point: (%g,%g)\n",value,X(ind,:));

%(8) Fill the feasible region and add labels
x=X(:,1);
y=X(:,2);
scatter(x,y,'*');
k=convhull(x,y);
fill(x(k),y(k),'m');
hold on
xlim([0 max(x)+1]);
ylim([0 max(y)+1]);

%(9) Plot to check
x1=0:0.1:max(B);
for z=1:8:value
    y=(z-c(1)*x1)/c(2);
    plot(x1,y);
    hold on
    drawnow;
    pause(0.001);
end
hold on

%% STANDARD LPP
%% Convert the problem into standard form 
% Max z = 2x1-3x2+6x3 
% s.t.x1-3x3>=4 
%     2x1-8x2+3x3<=4 
%     x1+x2>=-7 
%     x1,x2,x3>=0 

clc
clear all
format short

%(1) Define the function
c=[2,-3,6];
A=[1,0,-3;2,-8,3;1,1,0];
B=[4;4;-7];
sign=[-1;1;-1];
n=size(A,1);
s=eye(n);

%(2) Change the constraints according to the RHS
for i=1:n
    if B(i)<0
        A(i,:)=-A(i,:);
        B(i)=-B(i);
    end
end

%(3) Change the s values for the signs
for i=1:n
    if sign(i)<0
        s(i,:)=-s(i,:);
    end
end

%(4) display everything in form of table
obj_fun=array2table(c);
obj_fun.Properties.VariableNames(1:size(obj_fun,2))={'x1','x2','x3'}
Mat=[A s B]
std_lpp=array2table(Mat);
std_lpp.Properties.VariableNames(1:size(std_lpp,2))={'x1','x2','x3','s1','s2','s3','RHS'}

%% ALGEBRAIC METHOD
%% QUESTION 1  
% TO OBTAIN BFS USING ALGEBRAIC METHOD 
% Question 
% Max Z= 2x1+3x2+4x3+7x4  
% st: 2x1+3x2-x3+4x4=8 
% x1-2x2+6x3-7x4=-3 
% %xi>=0; i=1,2,3,4 

clc
clear all
format short
c=[2,3,4,7];
A=[2,3,-1,4;-1,2,-6,7];
B=[8;3];
objective=1;
[m,n]=size(A);
nab=nchoosek(n,m);
t=nchoosek(1:n,m);
sol=[];
if n>=m
    for i=1:nab
        y=zeros(n,1);
        X=(A(:,t(i,:)))\B;
        if all(X>=0&X~=inf&X~=-inf)
            y(t(i,:))=X;
            sol=[sol y];
        end
    end
    disp("Solution:");
    disp(sol)
else
    error("Number of variables is less than the number of constraints");
end

if any(X==0)
    disp("DEGENERATE");
else
    disp("NON-DEGENERATE");
end

Z=c*sol;
if(objective==1)
    [Zval,Zindex]=max(Z);
else
    [Zval,Zindex]=min(Z);
end
BFS=sol(:,Zindex);
[Optimal_value]=[BFS' Zval];
Optimal_bfs=array2table(Optimal_value);
Optimal_bfs.Properties.VariableNames(1:size(Optimal_bfs,2))={'x1','x2','x3','x4','val'};
disp(Optimal_bfs)

%% REPEAT
% TO OBTAIN BFS USING ALGEBRAIC METHOD 
% Question 
% Max Z= 2x1+3x2+4x3+7x4  
% st: 2x1+3x2-x3+4x4=8 
% x1-2x2+6x3-7x4=-3 
% %xi>=0; i=1,2,3,4 

c=[2,3,4,7];
A=[2,3,-1,4;-1,2,-6,7];
B=[8;3];
objective=1;
[m,n]=size(A);
nab=nchoosek(n,m);
t=nchoosek(1:n,m);
sol=[];
if(n>=m)
    for i=1:nab
        y=zeros(n,1);
        X=A(:,t(i,:))\B;
        if(X>=0&X~=inf&X~=-inf)
            y(t(i,:))=X;
            sol=[sol y];
        end
    end
else
    error("Number of variables is less than the number of constraints");
end
disp("Solution:");
disp(sol);

if any(X==0)
    disp("DEGENERATE");
else
    disp("NON-DEGENERATE");
end

Z=c*sol;
if(objective==1)
    [Zval,Zindex]=max(Z);
else
    [Zval,Zindex]=min(Z);
end
BFS=sol(:,Zindex);
[Optimal_solution]=[BFS' Zval];
opt=array2table(Optimal_solution);
opt.Properties.VariableNames(1:size(opt,2))={'x1','x2','x3','x4','bfs_val'};
disp(opt)
        
%% TEST CODE
%% Convert the problem into standard form 
% Max z = 2x1-3x2+6x3 
% s.t.x1-3x3>=4 
%     2x1-8x2+3x3<=4 
%     x1+x2>=-7 
%     x1,x2,x3>=0 

clc
clear all
format short

%(1) Define the function
c=[3,2,-1,1];
A=[1 2 1 -1;-2 -4 1 1];
B=[5;-1];
sign=[-1;1;-1];
n=size(A,1);
s=eye(n);

%(2) Change the constraints according to the RHS
for i=1:n
    if B(i)<0
        A(i,:)=-A(i,:);
        B(i)=-B(i);
    end
end

%(3) Change the s values for the signs
for i=1:n
    if sign(i)<0
        s(i,:)=-s(i,:);
    end
end

%(4) display everything in form of table
obj_fun=array2table(c);
obj_fun.Properties.VariableNames(1:size(obj_fun,2))={'x1','x2','x3'}
Mat=[A s B]
std_lpp=array2table(Mat);
std_lpp.Properties.VariableNames(1:size(std_lpp,2))={'x1','x2','x3','s1','s2','s3','RHS'}

