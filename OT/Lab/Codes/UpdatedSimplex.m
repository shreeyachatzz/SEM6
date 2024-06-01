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