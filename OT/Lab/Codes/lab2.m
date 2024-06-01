% %% Point of intersection with axes and between lines
% clc
% clear all
% format short
% %Max Z=x1+2x2
% %s.t. x1+x2<=5
% %     2x1+3x2<=10
% %     x1,x2>=0
% %     x1,x2=0, to find intersection points
% %Input parameters
% A= [1,1;2,3;0,1;1,0];
% B=[5;10;0;0];
% C=[1,2]; %cost of objective function
% pt=[0;0]; %initial start
% size_row= size(A,1);
% size_column= size(A,2);
% 
% %A1=[1,1], B1=[5]; A1X=B1
% %A2=[2,3], B2=[10]; A2X=B2
% %[A1;A2]X=[B1;B2]
% %A3X=B3
% for i=1:size_row
%     A1=A(i,:);
%     B1=B(i,:); %B(i)= B(i,:); since B is column matrix
%     for j=i+1:size_row
%         A2=A(j,:);
%         B2=B(j,:); %B(j)= B(j,:); since B is column matrix
%         A3=[A1;A2];
%         B3=[B1;B2];
%         X=inv(A3)*B3;
%         pt=[pt  X];
%     end
% end
% 
% X=pt'; %tranpsose for better view
% X=unique(X,'rows') %to make sure there are no repetition


%% QUESTION 2: Graphical method to solve 
% Max Z= 2x1+x2
% x1+2x2<=10
% x1+x2<=6
% x1-2x2<=1
% x2=0
% x1=0
clc
clear all
format short
%INPUT PARAMETERS
c=[5,4]; %cost objective function
A=[1,2;1,1;1,-2;0,1;1,0]; 
B=[10;6;1;0;0];
n=size(A,1);
x1=0:0.01:max(B);

for i=1:n-2 %we take n-2 since we are also taking x1=0 and x2=0 as they have no significance in our graph
    y(i,:)=(B(i)-A(i,1)*x1)/A(i,2);
end

%DRAWING THE LINES
for i=1:n-2
    y(i,:)=max(0,y(i,:));
    plot(x1,y(i,:),'linewidth',4)
    hold on
end
hold on
%FINDING THE POINT OF INTERSECTION
pt=[0;0];
for i=1:size(A,1)
    A1=A(i,:);
    B1=B(i,:);
    for j=i+1:size(A,1)
        A2=A(j,:);
        B2=B(j,:);
        A3=[A1;A2];
        B3=[B1;B2];
        %X3=inv(A3)*B3
        X3=A3\B3;
        if(X3>=0)%since the number of chairs can never be negative
            pt= [pt X3];
        end
    end
end

X=pt';
X=unique(X,'rows')%solution
hold on

% KEEP ONLY FEASIBLE POINTS
x1=X(:,1);
x2=X(:,2);

for i=1:n-2
    ind=find(A(i,:)*X'>B(i))
    X(ind,:)=[]
end

% EVALUATE THE OBJECTIVE FUNCTION VALUE
obj_val=c*X';
[value, ind]=max(obj_val);
value;
X(ind,:)
Optimal=[X(ind,:) value];

%% Shaded feasible region

x=X(:,1);
y=X(:,2);
scatter(X(:,1),X(:,2),'*')
hold on
k=convhull(x,y)%the shaded region where a and y is satisfied
fill(x(k),y(k),'m')

%% setting the axes
xlim([0 max(x)+1])
ylim([0 max(y)+1])

xlabel('x-axis')
ylabel('y-axis')
title('Feasible region of the linear programming problem')
legend('x_1+2x_2\leq10')
