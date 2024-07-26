%   Bayesian nonnegative CP factorization for tensor completion
%   mode inference (learning) via Gibbs sampling 
% 
%   coding by Zecan Yang in Aug. 2022
%   E-mail: zecanyang@gmail.com

% clc
close all;
clear all;

% % add code path
addpath(genpath('./lib'));

%% Generating synthetic data
I = 30; 
J = 40;
K = 50;
r = 5;
MaxRank = r + 15;
randn('state',1); rand('state',1); %#ok<RAND>
orgTensor = double(ktensor({rand(I, r), rand(J, r), rand(K, r)}));


sigma = 1e-3; % set noise 
miss_rate = 0.5;

% get binary tensor, if==1, observed; else, is missing
Is = size(orgTensor);
observe_idx = sparse_miss(Is, miss_rate, 'rand');

% generae Noise tensor;
NoiseT = sigma*randn(Is); 

% generate observe data;
yTensor = orgTensor.*observe_idx + NoiseT.*(1-observe_idx);


% get observation index
Map = observe_idx; % binary tensor

%% set hyper-parameters for model;
hyperparameters.alpha_tau_0 = 1e-6;
hyperparameters.beta_tau_0 = 1e-6;
hyperparameters.alpha_lambda_0 = 1e-6;
hyperparameters.beta_lambda_0 = 1e-6;


hyperparameters.orgTensor = orgTensor;
hyperparameters.NoiseT = NoiseT;
hyperparameters.DIMRED = true;
ARD = true;
train_iter = 30;
burn_out_iter = 10;

%% train model
bntf_gibbs = BayesianNNCP_Gibbs(yTensor, Map, MaxRank, ARD, hyperparameters); % get model
bntf_gibbs = bntf_gibbs.initialize();  % init model
bntf_gibbs = bntf_gibbs.run(train_iter, burn_out_iter);  % train model

