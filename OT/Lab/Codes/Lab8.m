%% GRADIENT DESCENT METHOD
%% QUESTION
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

