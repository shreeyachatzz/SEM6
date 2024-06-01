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