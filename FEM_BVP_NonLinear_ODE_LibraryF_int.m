%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FEM Code for Non-linear 1D BVP %%%%%%%%%
% u''+uu'-u=exp(2x), u(0)=1, u(1)=e %%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function FEM_BVP_NonLinear_ODE_LibraryF_int()
clear all; close all; clc
n=5;                        % No. of Element
nn=n+1;                     % No. of Nodes
lgth=1;                     % Domain length
he=lgth/n;                  % Length of each element
x=[0:he:lgth];              % Data points for independent variable
AC=0.00005;                 % Accuracy
F=zeros(nn,1);              % Initialization
F(1)=exp(0); F(nn)=exp(1);  % Boundary Conditions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Direct Iteration Process to handle nonlinear problem %%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c=1.0;
count=0;                    % Initialization for count for iterations
tic                         % Time start
while(c>0)
    [F1]=assembly(F,n,he);
    c=0.0;
    for i=1:nn
        if(abs(F(i)-F1(i))>AC)
            c=c+1;
            break;
        end
    end
    F=F1;
    count=count+1;
end
disp('Hence Solution');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output for primary and secondary variables %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
diff=abs(F-exp(x)');

fprintf('No of elements=%d\n',n)
disp('         x     FEM   Exact    Error')
disp([x',F,exp(x)',diff])
fprintf('No of iterations=%d\n',count)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ploting of primary variable %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot(x,F,'--rs','LineWidth',2)
xlabel('x')
ylabel('u(x)')
title('Solution Plot to Given BVP')

toc                                 % Give total time
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Derivation of element matrix and assembly %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [F1]=assembly(F,n,he)
nn=n+1;
K=zeros(nn,nn);                 % Initialization of main matrix 
R=zeros(nn,1);                  % Initialization of RHS matrix
syms x                          % x as symbolic variable
S=[1-x/he,x/he];                % Linear shape function
dS=diff(S,x);                   % Differentiation of shape function
lmm=[];
for i=1:n
    lmm=[lmm;[i,i+1]];          % Connectivity matrix
end

for i=1:n
    lm=lmm(i,:);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Generation of element matrix (k11) and RHS matrix (f1) %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    k11=-int(dS'*dS,x,0,he)+(int(S'*dS*S(1),x,0,he)*F(lm(1))...
        +int(S'*dS*S(2),x,0,he)*F(lm(2)))-int(S'*S,x,0,he);
    
    f1=int(exp(2*(x+(i-1)*he))*S',x,0,he);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Assembly according to connectivity matrix %%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    K(lm,lm)=K(lm,lm)+k11;
    R(lm)=R(lm)+f1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Imposition of boundary conditions %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
K(1,:)=0.0;   K(nn,:)=0.0;
K(1,1)=1.0;   K(nn,nn)=1.0;
R(1,1)=F(1);  R(nn,1)=F(nn);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solution of systems of equations (F1) %%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d=K\R;
F1=d;
end    