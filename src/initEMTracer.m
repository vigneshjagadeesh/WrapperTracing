%% Topology CRF ... Code written for manuscript submitted to TMI 
% Three inference techniques are implemented
% Greedy Inference - 1
% DP Inference - 2
% PF Inference - 3
% Author: Vignesh Jagadeesh
% Last Modified: 2/6/2011
clear; close all; clc;
addpath( genpath( '../' ) );
%% Preprocess the data
load splitMerge0.mat;
% Hyperparameters
intRange =  [50 100 150]'; %100
nlinkRange = [5 10 15]'; %30
inference = 1;
USEFLOW = 0;

% Initializations Common to all methods
[noRows noCols numSlices] = size(I_stack);
[edges, no_nlinks] = vrl_Construct2DLattice([noRows noCols],1);
% [edges no_nlinks] = vrlConstructGraph([noRows noCols], 1);
currOptima = zeros( size(I_stack) );
max_dist = sqrt(noRows^2 + noCols^2);    
numParams = numel(nlinkRange) * numel(intRange);
currProposals = cell( numParams, 1 );
thetaOpt = [];
segFgPrior = roipoly( I_stack(:,:,1) ); 
save('userInit.mat', 'segFgPrior');
%load('userInit.mat');

currOptima(:,:,1) = segFgPrior;
intShape = repmat(intRange', numel(nlinkRange), 1);
theta = [repmat(nlinkRange,numel(intRange),1) intShape(:)];


% Method Specific Initialization
if(inference == 1)
    qualityMat = zeros(numParams, 1); 
elseif(inference == 2)
    delDP = zeros(numParams, numSlices);
    delDP(:, 1) = 1/numParams;
    phiDP = zeros(numParams, numSlices);
    obsDP = ones( numParams, numSlices );
    qualityMatrix = zeros(numParams, numParams);
    for iter = 1:numParams
        prevProposals{iter} = segFgPrior; %#ok<SAGROW>
    end
elseif(inference == 3)
    numParams = 4;
    qualityMat = zeros(numParams, 1); 
    theta = theta( ceil(linspace(1,size(theta,1),numParams) ), : );
    wPF = ones( numParams, 1 ) / 100; % Generate uniform distribution
    tvPF = [3 3];
end