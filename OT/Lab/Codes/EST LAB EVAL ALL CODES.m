%% EST LAB EVAL
%% SIMPLEX METHOD
% Max Z= -x1+3x2-2x3
% s.t.  3x1-x2+2x3<=7
%       -2x1+4x2<=12
%       -4x1+3x2+8x3<=10
%       xi>=0  i=1-3
clc
clear all
format short
% Input Phase
Variables = {'x1','x2','x3','s1','s2','s3','sol'};
Cost = [-1 3 -2 0 0 0 0];
a=[3 -1 2; -2 4 0; -4 3 8];
s=eye(size(a,1));
b=[7;12;10];
A=[a s b];
s=eye(size(A,1));

%FINDING STARTING BFS
BV=[];
for j=1:size(s,2)
    for i=1:size(A,2)
        if A(:,i)==s(:,j)
            BV=[BV i];
        end
    end
end

% COMPUTE VALUE OF TABLE
B= A(:,BV);
A= inv(B)*A;
ZjCj= Cost(BV)*A-Cost;

% TO PRINT THE TABLE
fprintf('Simplex Table to solve: \n')
ZCj = [ZjCj;A];
SimpTable = array2table(ZCj);
SimpTable.Properties.VariableNames(1:size(ZCj,2))=Variables;
disp(SimpTable);

% SIMPLEX METHOD START
RUN =true;
while RUN
    ZC = ZjCj(:,1:end-1);
    if any(ZC<0)
        fprintf('Current BFS is NOT OPTIMAL\n');
        [Entval,pvt_col]=min(ZC);
        fprintf('Entering Column = %d \n',pvt_col);
        %finding leaving var
        sol = A(:,end);
        Column = A(:,pvt_col);
        if all(Column)<=0
            fprintf('Solution is UNBOUNDED');
        else
            for i=1:size(Column,1)
                if Column(i)>0
                    ratio(i)=sol(i)./Column(i);
                else
                    ratio(i)=inf;
                end
            end
            [minR, pvt_row]=min(ratio);
            fprintf('Leaving Row = %d\n',pvt_row);
            % UPDATE THE BV & TABLE
            BV(pvt_row)=pvt_col;
            % Updating the table
            B=A(:,BV);
            A= inv(B)*A;
            ZjCj = Cost(BV)*A-Cost;
            %to print intermediate table
            fprintf('Table after iteration: \n')
            ZCj = [ZjCj;A];
            TABLE = array2table(ZCj);
            TABLE.Properties.VariableNames(1:size(ZCj,2))=Variables;
            disp(TABLE);
        end
    else
        RUN = false;
        fprintf('CURRENT BFS IS OPTIMAL \n');
    end
end

%FINAL OPTIMAL SOLUTION PRINT:
% TO PRINT THE TABLE
FINAL_BFS= zeros(1,size(A,2));
FINAL_BFS(BV) = A(:,end);
FINAL_BFS(end) = sum(FINAL_BFS.*Cost);

% TO PRINT THE TABLE
OptimalBFS = array2table(FINAL_BFS);
OptimalBFS.Properties.VariableNames(1:size(OptimalBFS,2))=Variables;
fprintf("Final Optimal Table:\n");
disp(OptimalBFS);
val=(FINAL_BFS(end));
fprintf("Optimal Value of Z: %0.2f\n",val);

%% BIG-M METHOD
%% QUESTION 1
% Solve LPP using simplex using Simplex Algorithm with Big-M method
% Minimize Z = 2x1+x2
% s.t.  3x1 + x2 = 3
%       4x1 + 3x2 >= 6
%       x1 + 2x2 <= 3
%       xi >= 0  i=1-3
% Maximize Z = -2x1-x2-Ma1-Ma2
% s.t.  3x1 + x2 +a1 = 3
%       4x1 + 3x2 -s2 + a2 = 6
%       x1 + 2x2 +s3 = 3
%       xi >= 0  i=1-3
clc
clear all
format short
% Input Phase
Variables = {'x1','x2','s2','s3','a1','a2','Sol'};
M=1000;
Cost = [-2,-1,0,0,-M,-M,0];
a=[3,1,0,0,1,0; 4,3,-1,0,0,1;1,2,0,1,0,0];
b=[3;6;3];
A=[a b];
s=eye(size(A,1));

%FINDING STARTING BFS
BV=[];
for j=1:size(s,2)
    for i=1:size(A,2)
        if A(:,i)==s(:,j)
            BV=[BV i];
        end
    end
end

% COMPUTE VALUE OF TABLE
B= A(:,BV);
A= inv(B)*A;
ZjCj= Cost(BV)*A-Cost;

% TO PRINT THE TABLE
fprintf('Simplex Table to solve: \n')
ZCj = [ZjCj;A];
SimpTable = array2table(ZCj);
SimpTable.Properties.VariableNames(1:size(ZCj,2))=Variables;
disp(SimpTable);

% SIMPLEX METHOD START
RUN =true;
while RUN
    ZC = ZjCj(:,1:end-1);
    if any(ZC<0)
        fprintf('Current BFS is NOT OPTIMAL\n');
        [Entval,pvt_col]=min(ZC);
        fprintf('Entering Column = %d \n',pvt_col);
        %finding leaving var
        sol = A(:,end);
        Column = A(:,pvt_col);
        if all(Column)<=0
            fprintf('Solution is UNBOUNDED');
        else
            for i=1:size(Column,1)
                if Column(i)>0
                    ratio(i)=sol(i)./Column(i);
                else
                    ratio(i)=inf;
                end
            end
            [minR, pvt_row]=min(ratio);
            fprintf('Leaving Row = %d\n',pvt_row);
            % UPDATE THE BV & TABLE
            BV(pvt_row)=pvt_col;
            % Updating the table
            B=A(:,BV);
            A= inv(B)*A;
            ZjCj = Cost(BV)*A-Cost;
            %to print intermediate table
            fprintf('Table after iteration: \n')
            ZCj = [ZjCj;A];
            TABLE = array2table(ZCj);
            TABLE.Properties.VariableNames(1:size(ZCj,2))=Variables;
            disp(TABLE);
        end
    else
        RUN = false;
        fprintf('CURRENT BFS IS OPTIMAL \n');
    end
end

%FINAL OPTIMAL SOLUTION PRINT:
% TO PRINT THE TABLE
FINAL_BFS= zeros(1,size(A,2));
FINAL_BFS(BV) = A(:,end);
FINAL_BFS(end) = sum(FINAL_BFS.*Cost);

% TO PRINT THE TABLE
OptimalBFS = array2table(FINAL_BFS);
OptimalBFS.Properties.VariableNames(1:size(OptimalBFS,2))=Variables;
fprintf("Final Optimal Table:\n");

disp(OptimalBFS);
val=(FINAL_BFS(end));
fprintf("Optimal Value of Z: %0.2f\n",val);

%% LEAST COST METHOD
clc
clear all
format short
% Input the parameters
Cost=[11,20,7,8;21,16,10,12;8,12,18,19];%Cost matrix
S=[50,40,70];%Supply
D=[30,25,35,40];%Demand

% To check if UBTP and fix it if UBTP
if(sum(S)==sum(D))
    fprintf("This is a Balanced Transportation Problem\n")
else
    fprintf("This is an Unbalanced Transportation Problem\n")
    if(sum(S)<sum(D))
        % Add a dummy column with size matching the number of origins (in Demand)
        Cost(:,end+1) = zeros(1,size(D,2));
        % Adding an element in S
        S(end+1)=sum(D)-sum(S);
    else
        % Add a dummy row with size matching the number of destinations (in Supply)
        Cost(:,end+1) = zeros(1,size(S,2));
        %Adding an element in D
        D(end+1)=sum(S)-sum(D);
    end
end

% START LCEM to find MINIMUM COST and MAXIMUM ALLOCATION
ICost=Cost;%make a copy of Cost to make edits during the process
X= zeros(size(Cost));%initial allocation to all places = 0
[m,n]=size(Cost);% to find the number of rows and columns to calcumate m+n-1
nBFS=m+n-1; %total number of BFS cells
% Now find cells with min cost and max allocation
for i=1:size(Cost,1)
    for j=1:size(Cost,2)
        hh=min(Cost(:)); %to find the min cost value in the entire matrix
        [row_ind,col_ind]=find(hh==Cost); %to find the index where hh is equal to the value in cost(pos of minimum cell)
        x11=min(S(row_ind),D(col_ind));%finding the minimum values for allocation for the selected col/row index
        [val,ind]=max(x11);%Choose the max among the min values for maximum allocation and also find it's index
        ii = row_ind(ind); %to identify the row position in the Cost matrix
        jj = col_ind(ind); %to identify the column poisiton in the Cost matrix
        y11=min(S(ii),D(jj)); %to find the value to be inserted into X
        X(ii,jj)=y11; %allocating the value at X at (ii,jj)
        Cost(ii,jj)=Inf; %Allocate a maximum value at the cell where we just found the BV so that it is not considered in further iterations
        %updating the S and D values
        S(ii)=S(ii)-y11;
        D(jj)=D(jj)-y11;
    end
end

% To print the initial BFS
fprintf("Initial BFS:\n")
IBFS=array2table(X);
IBFS.Properties.VariableNames(1:size(IBFS,2))={'D1','D2','D3','D4','D5'};
disp(IBFS)

%%TO CHECK DEGENERATE
FinalBFS=length(nonzeros(X));
if(FinalBFS==nBFS)
    fprintf('The initial BFS is NOT DEGENERATE\n');
else
    fprintf('The initial BFS is DEGENERATE\n');
end

%%TO COMPUTE IBFS TP Cost
TPCost=sum(sum(ICost.*X));
fprintf('The initial Transportation Cost is = %d\n', TPCost)

%% GRADIENT DESCENT METHOD
% QUESTION
% Min f(x1,x2)=x1-x2+2x1^2+2x1x2+x^2;
% Starting at x1=(1,1)
clc
clear all
format short

% PHASE1- INPUT THE PARAMETER
syms x1 x2
f1=x1-x2+2*x1^2+2*x1*x2+x2^2;
fx=inline(f1); %to convert the text into a function
fobj=@(x) fx(x(:,1), x(:,2)); %to find the value of the function at given. x(:,1)-> picks the first column of the innput(x1,x2) and x(:,2) picks the second column

% PHASE2 - TO COMPUTE THE GRADIENT
grad=gradient(f1);
G=inline(grad);
Gradx=@(x) G(x(:,1), x(:,2));

%PHASE 3- COMPUTE HESSIAN
H1=hessian(f1);
Hx=inline(H1);% the hessian matrix is converted into a inline matrix(reshaping is done)

%PHASE 4: ITERATION
x0=[1,1]; %initial value (for this iteration)
max_itr=4;%(this maximum number of iteration will be given in the question or some other terminating condition will be given)
tol=10^(-3);%error tolerance(given)
iter=0; %initial counter/iterator
%X=[]; %initial empty vector where we will be storing the values

%norm(x1,x2)=sqrt(x1^2+x2^2)
while (norm(Gradx(x0))>tol && iter<max_itr)
%       X=[X x0]; %initialize this vector with the initial pointer
        S=-Gradx(x0); %compute gradient
        H=Hx(x0); %compute hessian at this point
        lambda=(S'*S)./(S'*H*S); %computing the value of lambda
        xnew=x0+lambda.*S'; %updated value after new iteration
        x0=xnew; %save new x after this iteration
        iter=iter+1; %update iterator
end

% PHASE 5 PRINT THE SOLUTION
fprintf('Optimal Solution x=[%f , %f]\n',x0(1),x0(2))
fprintf('Optimal Value f(x)=%f \n',fobj(x0))