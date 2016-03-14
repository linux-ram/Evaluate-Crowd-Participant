function prediction=MLLabeler_testing(Xtst,alpha)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function prediction=MLLabeler_testing(Xtst,alpha)
% Y.Yan Oct. 20th 2010    
% this function is the testing file you should call
% only after param alpha is acquired
% note that testing data Xtst should follow the same 
% dimensions and scales as the training data you used
% when training alpha
% INPUT:
% Xtst: testing data
% alpha: classifier acquired by calling function MLLabeler 
%
% OUPUT
% prediction: probabilities of binary hidden true labels of testing data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

prediction=1./(1+exp(-1*([Xtst ones(size(Xtst,1),1)]*alpha)));

return