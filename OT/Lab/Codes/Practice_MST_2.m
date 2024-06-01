%% TUT 2: 1(i)
% Max Z= 5x1+3x2
% st
% 3x1+5x2<=15
% 5x1+2x2<=10
% x1,x2>=0
clf
clc
clear all
format short
%(1) Declaring the function
c=[5,3];
A=[3,5;5,2;0,1;1,0];
B=[15;10;0;0];
const=[1,1];
objective=1;
n=size(A,1);
x1=0:0.01:max(B);

%(2) Finding the values of y
for i=1:n-2
    y(i,:)=(B(i)-A(i,1)*x1)/A(i,2);
end

%(3) Plotting the lines
for i=1:n-2
    y(i,:)=max(0,y(i,:));
    plot(x1,y(i,:), "LineWidth",4);
    hold on;
end
hold on;
%%
%(4) Finding the points of intersection
pt=[0;0];
for i=1:n
    A1=A(i,:);
    B1=B(i);
    for j=i+1:n
        A2=A(j,:);
        B2=B(j);
        A3=[A1;A2];
        B3=[B1;B2];
        X3=inv(A3)*B3;
        if X3>=0
            pt=[pt X3];
        end
    end
end
%(5) inverse pt and store unique points only
X=pt';
X=unique(X,'rows')
hold on 

%(6) Finding the feasible solutions
for i=1:n-2
    if (const(i)==1)
        ind=find(A(i,:)*X'>B(i));
        X(ind,:)=[];
    else
        ind=find(A(i,:)*X'<B(i));
        X(ind,:)=[];
    end
end

%(7) Finding optimal solution
obj_val=c*X';
if(objective==1)
    [value,ind]=max(obj_val);
else
    [value,ind]=min(obj_val);
end
fprintf("Optimal Value: %f\n",value);
fprintf("Optimal Point: %g %g\n",X(ind,:));

%(8) Shading feasible region
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
for z=0:0.2:value
    y=(z-c(1)*x)/c(2);
    plot(x,y);
    hold on
    drawnow
    pause(0.001)
end
hold on

%% CONVERT TO STANDARD FORM

c=[2,3];
A=[1,3;3,2];
B=[6;6];
sign=[1;1];
[n,m]=size(A);
s=eye(n);

for i=1:n
    if B(i)<0
        A(i,:)=-A(i,:);
        B(i)=-B(i);
        sign(i)=-sign(i);
    end
end

for i=1:n
    if sign(i)<0
        s(i,:)=-s(i,:);
    end
end

obj_fun=array2table(c);
obj_fun.Properties.VariableNames(1:size(obj_fun,2))={'x1','x2'}
Mat=[A s B]
std=array2table(Mat);
std.Properties.VariableNames(1:size(std,2))={'x1','x2','s1','s2','rhs'}

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
[m,n]=size(A);
nab=nchoosek(n,m);
t=nchoosek(1:n,m);
objective=1;

sol=[];
if n>=m
    for i=1:nab
        y=zeros(n,1);
        val=(A(:,t(i,:)))\B;
        if (val>=0 & val~=-inf & val~=inf)
            y(t(i,:))=val;
            sol=[sol y];
        end
    end
else
    error("Number of variables less than number of constrains");
end

disp("Solution:");
disp(sol);

if any(val==0)
    disp("DEGENERATE");
else
    disp("NON-DEGENERATE");
end

obj=c*sol;
if objective==1
    [Zval,Zind]=max(obj);
else
    [Zval,Zind]=min(obj);
end

BFS=sol(:,Zind);
[Optimal_solution]=[BFS' Zval];
opt=array2table(Optimal_solution);
opt.Properties.VariableNames(1:size(opt,2))={'x1','x2','x3','x4','bfs_val'};
disp(opt)

