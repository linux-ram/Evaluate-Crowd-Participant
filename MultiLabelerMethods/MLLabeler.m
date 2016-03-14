function [alpha W PtildZ Pz_x maxLogScr Py_xz Py_xzbar]=MLLabeler(X, Y, MChoice, alpha_ini, W_ini)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [alpha W PtildZ Pz_x maxLogScr]=MLLabeler(X, Y, MChoice, alpha_ini, W_ini)
% Y.Yan Oct. 20th 2010        
% INPUT
% MChoice: 1 Gaussian; 2 Binomial
% yi should be taking binary values of {0,1}
% X should be in the form of N*D (N:no of samples; D: feature dimensions)
% we make the alpha form as Xext*alpha
% W should be the model parameter of each annotator
% alpha_ini and W_ini are optional params which are for the 
% initialization of classifier alpha and performance param W
% specifically when alpha_ini=1, it means to randomly initialize alpha
% instead of setting alpha equals to 1 directly.
%
% OUTPUT
% alpha: the classifier you probably care the most
% W: performance measures for each labeler
% PtildZ: intermediate params used when processing EM
% Pz_x: true label prediction for training data
% maxLogScr: maximum loglikelihood score after processing EM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


N=size(X,1);
D=size(X,2);
T=size(Y,2);
Xext=[X ones(N,1)];
epsilon=1e-8;
Zcrs=max(sum(Y,2)./T,1e-4*ones(N,1));

if nargin<3
    MChoice=1; % gaussian
end

if (nargin<4)||(isempty(alpha_ini))
    % initial alpha
    alpha_ini=(Xext'*Xext+norm(Xext)*1e-4*eye(size(Xext,2),size(Xext,2)))\(Xext'*(-1.*(log(1./Zcrs-1+1e-4))));
end

if (nargin>=4)&&(size(alpha_ini,1)==1)
    if alpha_ini==1;
        clear alpha_ini
        alpha_ini=randn(D+1,1);
    end
end

if nargin<5
    % initialize W
    W_ini=zeros(D+1,T);
end

%%

alpha=zeros(D+1,1);
alpha_new=alpha_ini;

W=zeros(D+1,T);
W_new=W_ini;

% PtildZ initialize
PtildZ=max(sum(Y,2)/T, 1e-5);

% E-M algorithm
disp('E_M algorithm for alpha and W starts now:');
step=1;
while(norm(alpha_new-alpha)^2+norm(W-W_new, 'fro')^2>=epsilon)
    clc
    disp(strcat('step: ', num2str(step)));
    pause(.2);
    
    alpha=alpha_new;
    W=W_new;
    
    % M step
    % update alpha
    disp('updating alpha...');
    maxLogScr=0;
    %...function from M.Schmidt
    options.Method = 'lbfgs';
    %options.Method = 'cg';
    options.Display = 'none';
    [alpha_new maxLogScrPlus]=minFunc(@(alpha_func) MLfunction_alpha(X, Y, alpha_func, PtildZ), alpha, options);
    %--
    %...matlab toolbox
    %options=optimset('GradObj','on','Display','off', 'HessUpdate', 'steepdesc');
    %[alpha_new maxLogScrPlus]=fminunc(@(alpha_func) MLfunction_alpha(X, Y, alpha_func, PtildZ), alpha, options);
    %--
    
    maxLogScr=maxLogScr+maxLogScrPlus;
    % update W columnwise
    disp('updating W by each user...');
    for t=1:1:T
        %...function from M.Schmidt
        [W_new(:,t) maxLogScrPlus]=minFunc(@(Wt_func) MLfunction_W(X, Y, Wt_func, t, PtildZ, MChoice), W(:,t), options);
        %--
        %...matlab toolbox
        %[W_new(:,t) maxLogScrPlus]=fminunc(@(Wt_func) MLfunction_W(X, Y, Wt_func, t, PtildZ, MChoice), W(:,t), options);
        %--
 
        maxLogScr=maxLogScr+maxLogScrPlus;
    end
    
    % E-M maximum likelihood computation
    disp(strcat('the total likelihood score is: ', num2str(-1*maxLogScr)));
    pause(1);
    
    % E step
    Pz_x=1./(ones(N,1)+exp(-1*Xext*alpha_new));
    Py_xz=zeros(N,T);
    Py_xzbar=zeros(N,T);
    for i=1:1:N
        for t=1:1:T
            if MChoice==1 % gaussian
                sigma_i_t=1/(1+exp(Xext(i,:)*W_new(:,t)));
                Py_xz(i,t)=normpdf(Y(i,t),1,sigma_i_t);
                Py_xzbar(i,t)=normpdf(Y(i,t),0,sigma_i_t);
            else % binomial
                p_i_t=1/(1+exp(Xext(i,:)*W_new(:,t)));
                Py_xz(i,t)=(1-p_i_t)^(1-Y(i,t))*p_i_t^Y(i,t);
                Py_xzbar(i,t)=(1-p_i_t)^Y(i,t)*p_i_t^(1-Y(i,t));
            end
        end
    end
    PtildZ=zeros(N,1);
    for i=1:1:N
        PtildZ(i,1)=1/(1+(prod((Py_xzbar(i,:)+1e-10)./(Py_xz(i,:)+1e-10))*(1-Pz_x(i,1))/Pz_x(i,1)));
    end

    
    step=step+1;
end

return

%% ML functions
function  [F_alpha Grad_alpha]=MLfunction_alpha(X, Y, alpha, PtildZ) 

N=size(X, 1);
D=size(X, 2);
T=size(Y, 2);
Xext=[X ones(N,1)];

Grad_alpha=zeros(D+1, 1);
F_alpha=0;

for i=1:1:N
    for t=1:1:T
        F_alpha=F_alpha-PtildZ(i,1)*log(1+exp(-1*Xext(i,:)*alpha));
        F_alpha=F_alpha-(1-PtildZ(i,1))*(Xext(i,:)*alpha+log(1+exp(-1*Xext(i,:)*alpha)));
    end     
    %Grad_alpha=Grad_alpha+(PtildZ(i,1)-1/(1+exp(-1*Xext(i,:)*alpha)))*Xext(i,:)';
    Grad_alpha=Grad_alpha+T*(2*PtildZ(i,1)-1)*exp(-1*Xext(i,:)*alpha)/((1+exp(-1*Xext(i,:)*alpha))^2)*Xext(i,:)';
end

Grad_alpha=-1*Grad_alpha;
F_alpha=-1*F_alpha;
return

function [F_W Grad_Wt]=MLfunction_W(X, Y, Wt, Tind, PtildZ, MChoice)

N=size(X, 1);
D=size(X, 2);

Grad_Wt=zeros(D+1,1);
F_W=0;
Xext=[X ones(N,1)];
for i=1:1:N
        if MChoice==1 % gaussian
            % function value
            sigma=1/(1+exp(Xext(i,:)*Wt));
            F_W=F_W+PtildZ(i,1)*log(max(normpdf(Y(i,Tind), 1, sigma), 1e-5));
            F_W=F_W+(1-PtildZ(i,1))*log(max(normpdf(Y(i,Tind), 0, sigma), 1e-5));
            
            Grad_Wt=Grad_Wt+((Y(i,Tind)^2-2*PtildZ(i,1)*Y(i,Tind)+PtildZ(i,1))/(sigma^2)-1)*(sigma-1)*Xext(i,:)';
        
        else %binomial
            P_xt=1/(1+exp(Xext(i,:)*Wt));
            F_W=F_W+PtildZ(i,1)*log(max(P_xt^Y(i,Tind)*(1-P_xt)^(1-Y(i,Tind)), 1e-5));
            F_W=F_W+(1-PtildZ(i,1))*log(max(P_xt^(1-Y(i,Tind))*(1-P_xt)^Y(i,Tind), 1e-5));                       
            Grad_Wt=Grad_Wt+(-1)^Y(i,Tind)*(1-2*PtildZ(i,1))*P_xt*(P_xt-1)*Xext(i,:)';
        end
end
F_W=-1*F_W;
Grad_Wt=-1*Grad_Wt;
return