clc
clear all
format short

% 1. Initializing the variable
cost=[10,2,20,11;12,7,9,20;4,14,16,18];
S=[15,25,10];
D=[5,15,15,15];

% 2. Checking if Balanced or Unbalanced TP
if sum(S)==sum(D)
    disp ("BALANCED TRANSPORTATION PROBLEM");
else
    disp ("UNBALANCED TRANSPORTATION PROBLEM");
    if (sum(S)<sum(D))
        cost(end+1,:)=zeros(1,size(D,2));
        S(end+1)=sum(D)-sum(S);
    else
        cost(:,end+1)=zeros(size(S,1),1);
        D(end+1)=sum(S)-sum(D);
    end
end

% 3. Starting LCM
icost=cost;
X=zeros(size(cost));

while any(cost(:) < inf)
    min_val = min(cost(:));
    [minr, minc] = find(cost == min_val,1);

    if isempty(minr) || isempty(minc)
        break; % No more feasible assignments
    end

    X(minr,minc)=min(S(minr),D(minc));

    S(minr)=S(minr)-X(minr,minc);
    D(minc)=D(minc)-X(minr,minc);
    cost(minr,minc)=inf;
end

% 4. Display the result
fprintf("Allocation Matrix:\n");
table=array2table(X);
table.Properties.VariableNames(1:size(X,2))={'D1','D2','D3','D4'};
table.Properties.RowNames(1:size(X,1))={'O1','O2','O3'};
disp(table);

total_cost = sum(sum(X .* icost));
fprintf("Initial BFS= %f\n",total_cost);
