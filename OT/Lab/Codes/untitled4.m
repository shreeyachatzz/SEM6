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
c=[6,11];
A=[2,1;1,2;0,1;1,0];
B=[104;76;0;0];
const=[1;1];
n=size(A,1);
objective=1;
x1=0:0.01:max(B);

for i=1:n-2
    y(i,:)=(B(i)-(A(i,1)*x1))/A(i,2);
end

for i=1:n-2
    y(i,:)=max(0,y(i,:));
    plot(x1,y(i,:),"LineWidth",4);
    hold on
end
hold on

pt=[0;0];
for i=1:n
    A1=A(i,:);
    B1=B(i,:);
    for j=1:n
        A2=A(j,:);
        B2=B(j,:);
        A3=[A1;A2];
        B3=[B1;B2];
        X3=A3\B3;
        if X3>=0
            pt=[pt X3];
        end
    end
end

X=pt';
X=unique(X,'rows');

for i=1:n-2
    if const(i)==1
        ind=find(A(i,:)*X'>B(i));
        X(ind,:)=[];
    else
        ind=find(A(i,:)*X'<B(i));
        X(ind,:)=[];
    end
end

Z=c*X';

if(objective==1)
    [val,ind]=max(Z);
else
    [val,ind]=min(Z);
end

fprintf("Optimal Solution: %f\n",val);
fprintf("Optimal Point: %g %g\n",X(ind,:));

x=X(:,1);
y=X(:,2);
scatter(x,y,'yellow','*');
hold on;
k=convhull(x,y);
fill(x(k),y(k),'m');
xlim([0 max(x)+1]);
ylim([0 max(y)+1]);

x1=0:0.01:max(B);
for z=0:8:val
    y=(z-c(1)*x1)/c(2);
    plot(x1,y);
    hold on
    drawnow
    pause(0.01);
end
hold on

%% Convert the problem into standard form
% Max z = 3x1+5x2
% s.t. x1+2x2<=20
%     x1+x2<=15
%     x2>=6
%     x1,x2>=0
clc
clear all
format short

c=[3,5];
A=[1,2;1,1;0,1];
B=[20;15;6];
n=size(A,1);
s=eye(n);
sign=[1 1 -1];

for i=1:n
    if (B(i)<0)
        A(i,:)=-A(i,:);
        B(i)=-B(i);
        sign(i)=-sign(i);
    end
end

for i=1:n
    if(sign(i)<0)
        s(i,:)=-s(i,:);
    end
end

obj=array2table(c);
obj.Properties.VariableNames(1:size(obj,2))={'x1','x2'}
Mat=[A s B]
obj=array2table(Mat);
obj.Properties.VariableNames(1:size(obj,2))={'x1','x2','s1','s2','s3','rhs'}

%% Convert the problem into standard form 
% Max z = 2x1-3x2+6x3 
% s.t.x1-3x3>=4 
%     2x1-8x2+3x3<=4 
%     x1+x2>=-7 
%     x1,x2,x3>=0 

c=[2,-3,6];
A=[1,0,-3;2,-8,3;1,1,0];
B=[4;4;-7];
sign=[-1 1 -1];
n=size(A,1);
s=eye(n);

for i=1:n
    if(B(i)<0)
        B(i)=-B(i);
        s(i,:)=-s(i,:);
        A(i,:)=-A(i,:);
    end
end

for i=1:n
    if(sign(i)<0)
        s(i,:)=-s(i,:);
    end
end

obj=array2table(c);
obj.Properties.VariableNames(1:size(obj,2))={'x1','x2','x3'}
Mat=[A s B];
lpp=array2table(Mat);
lpp.Properties.VariableNames(1:size(lpp,2))={'x1','x2','x3','s1','s2','s3','rhs'}


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

c=[6,11];
A=[2,1;1,2;1,0;0,1];
B=[104;76;0;0];
sign=[1,1];
n=size(A,1);
x1=0:0.01:max(B);
objective=1;

for i=1:n-2
    y(i,:)=(B(i)-A(i,1)*x1)/A(i,2);
end

for i=1:n-2
    y(i,:)=max(0,y(i,:));
    plot(x1,y(i,:),'LineWidth',4);
    hold on
end
hold on

pt=[0;0];
for i=1:n
    A1=A(i,:);
    B1=B(i,:);
    for j=i+1:n
        A2=A(j,:);
        B2=B(j,:);
        A3=[A1;A2];
        B3=[B1;B2];
        X3=A3\B3;
        if X3>=0
            pt=[pt X3];
        end
    end
end

X=pt';
X=unique(X,'rows');

for i=1:n-2
    if(sign(i)==1)
        ind=find(A(i,:)*X'>B(i));
        X(ind,:)=[];
    else
        ind=find(A(i,:)*X'>B(i));
        X(ind,:)=[];
    end
end

Z=c*X';
if(objective==1)
    [val,ind]=max(Z);
else
    [val,ind]=min(Z);
end

fprintf("Optimal Value:%f\n",val);
fprintf("Optimal Solution: (%g,%g)\n",X(ind,:));

x=X(:,1);
y=X(:,2);
scatter(x,y,'yellow',"*");
k=convhull(x,y);
fill(x(k),y(k),'m');
xlim([0 max(x)+1]);
ylim([0 max(y)+1]);
hold on

x1=0:0.01:max(B);

for z=0:8:val
    y=(z-c(1)*x1)/c(2);
    plot(x1,y);
    hold on
    drawnow
    pause(0.01);
end

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

c=[2 3 4 7];
A=[2,3,-1,4;-1,2,-6,7];
B=[8;3];
[m,n]=size(A);
objective=1;

if n>m
    nab=nchoosek(n,m);
    t=nchoosek(1:n,m);
    
    sol=[];
    for i=1:nab
        y=zeros(n,1);
        X=(A(:,t(i,:)))\B;
        if(X>=0 & X~=-inf & X~=inf)
            y(t(i,:))=X;
            sol=[sol y];
        end
    end
else
    error ("Number of constraints more than the number of variables");
end

disp("Points:");
disp(sol);

val=c*sol;

if(objective==1)
    [val, ind]=max(val);
else
    [val,ind]=min(val);
end


Mat=[sol(:,ind)' val];
final=array2table(Mat);
final.Properties.VariableNames(1:size(final,2))={'x1','x2','x3','x4','value'}

if any(X==0)
    disp("DEGENERATE SOLUTION");
else
    disp("NON DEGENERATE SOLUTION");
end

%% QUESTION 3 
% TO OBTAIN BFS USING ALGEBRAIC METHOD
% Question
% Min Z= 5x2-2x1
% st: 2x1+5x2+s1=8
% x1+x2+s2=2
% %xi>=0
clc
clear all 
format short
c=[-2 5 0 0];
A=[2,5,1,0;1,1,0,1];
B=[8;2];
[m,n]=size(A);
objective=-1;

if(n>=m)
    nab=nchoosek(n,m);
    t=nchoosek(1:n,m);

    sol=[];
    for i=1:nab
        y=zeros(n,1);
        X=(A(:,t(i,:)))\B;
        if(X>=0 & X~=-inf & X~=inf)
            y(t(i,:))=X;
            sol=[sol y];
        end
    end
    disp("Points:");
    disp(sol);
else
    error("Number of constraints more than the number of variables");
end

Z=c*sol;

if(objective==1)
    [val, ind]=max(Z);
else
    [val, ind]=min(Z);
end

Mat=[sol(:,ind)' val];
final=array2table(Mat);
final.Properties.VariableNames(1:size(final,2))={'x1','x2','s1','s2','sol'}

%% QUESTION 2: Graphical method to solve 
% Max Z= 5x1+8x2
% 2x1+x2<=200
% x1+2x2<=150
% x2<=60
% x1,x2>=0

clc
clf
clear all
format short

c=[5,8];
A=[2,1;1,2;0,1;0,1;1,0];
B=[200;150;60;0;0];
const=[1;1;1];
objective=1;
n=size(A,1);
x1=0:0.01:max(B);

for i=1:n-2
    y(i,:)=(B(i)-A(i,1)*x1)/A(i,2);
end

for i=1:n-2
    y(i,:)=max(0, y(i,:));
    plot(x1,y,'LineWidth',4);
    hold on
end
hold on

pt=[0;0];
for i=1:n
    A1=A(i,:);
    B1=B(i,:);
    for j=1:n
        A2=A(j,:);
        B2=B(j,:);
        A3=[A1;A2];
        B3=[B1;B2];
        X3=A3\B3;
        if X3>=0
            pt=[pt X3];
        end
    end
end

X=pt';
X=unique(X,'rows');

for i=1:n-2
    if(const(i)==1)
        ind=find(A(i,:)*X'>B(i));
        X(ind,:)=[];
    else
        ind=find(A(i,:)*X'<B(i));
        X(ind,:)=[];
    end
end

Z=c*X';
if objective==1
    [val,ind]=max(Z);
else
    [val,ind]=min(Z);
end

disp("Optimal Solution")
disp(X(ind,:));
disp("Optimal Value:")
disp(val)


x=X(:,1);
y=X(:,2);
scatter(x,y,'yellow','*');
hold on;
k=convhull(x,y);
fill(x(k),y(k),'m');
hold on
xlim([0 max(x)+2]);
ylim([0 max(y)+2]);

for z=1:8:val
    y1=(z-c(1)*x1)/c(2);
    plot(x1,y1);
    hold on;
    drawnow;
    pause(0.001);
end

%% Convert the problem into standard form
% Min z = 40x1+24x2+12x3
% s.t. 20x1+50x2+10x3>=480
%     8x1+5x2+2x3<=72
%     4x1+5x2+3x3<=7
%     x1,x2,x3>=0

c=[40,24,12];
A=[20,50,10;8,5,2;4,5,3];
B=[480;72;7];
sign=[-1 1 1];
n=size(A,1);
s=eye(n);

for i=1:n
    if(B(i)<0)
        B(i)=-B(i);
        A(i,:)=-A(i,:);
        sign(i)=-sign(i);
    end
end

for i=1:n
    if sign(i)<0
        s(i,:)=-s(i,:);
    end
end

obj=array2table(c);
obj.Properties.VariableNames(1:size(obj,2))={'x1','x2','x3'}
Mat=[A s B];
lpp=array2table(Mat);
lpp.Properties.VariableNames(1:size(lpp,2))={'x1','x2','x3','s1','s2','s3','val'}


%% TO OBTAIN BFS USING ALGEBRAIC METHOD
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

sol=[];
if(n>=m)
    nab=nchoosek(n,m);
    t=nchoosek(1:n,m);

    for i=1:nab
        y=zeros(n,1);
        X=(A(:,t(i,:)))\B;
        if (X>=0 & X~=-inf & X~=inf)
            y(t(i,:))=X;
            sol=[sol y];
        end
    end
    disp("Points:")
    disp(sol)
else
    error("Number of varaibles less than number of constraints");
end

Z=c*sol;

if(objective==1)
    [val,ind]=max(Z);
else
    [val,ind]=min(Z);
end

Mat=[sol(:,ind)' val];
final=array2table(Mat);
final.Properties.VariableNames(1:size(final,2))={'x1','x2','x3','x4','val'};
final

%% TO OBTAIN BFS USING ALGEBRAIC METHOD
% Question
% Max Z= 2x1+3x2+4x3+7x4 
% st: 2x1+3x2-x3+4x4=8
% x1-2x2+6x3-7x4=-3
% %xi>=0; i=1,2,3,4

clc
clear all
format short

c=[2,3,4,7];
A=[2,3,-1,4;1,-2,6,-7];
B=[8;-3];
objective=1;
[m,n]=size(A);

sol=[];
if(n>=m)
    nab=nchoosek(n,m);
    t=nchoosek(1:n,m);
    for i=1:nab
        y=zeros(n,1);
        X=(A(:,t(i,:)))\B;
        if (X>=0 & X~=-inf & X~=inf)
            y(t(i,:))=X;
            sol=[sol y];
        end
    end
    disp("Points:");
    disp(sol);
else
    error("The number of variables is less than the number of constraints");
end

Z=c*sol;
if objective==1
    [val,ind]=max(Z);
else
    [val,ind]=min(Z);
end

mat=[sol(:,ind)' val];
p=array2table(mat);
p.Properties.VariableNames(1:size(p,2))={'x1','x2','x3','x4','val'}






    



         
