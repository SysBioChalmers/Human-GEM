function [] = generate_DepMap_models_ftINIT_Cluster(chunk, use1Plus1)
% generate tINIT models from RNA-Seq profiles
%
% Note: This function was designed specifically for use on the cluster
%
% Input:
%
%   chunk     integer ranging from 1 to 10, or 'test' to run a test that
%             generates only one model.
%
%   use1Plus1 if TRUE, run ftINIT with 1+1
%             (Default = FALSE)
%

if nargin < 2
    use1Plus1 = false;
end

% add paths
addpath(genpath('../components/RAVEN'));
addpath(genpath('../components/COBRA'));
addpath(genpath('../components/Human-GEM'));

cd '../components/Human-GEM/code/DepMapGeneEss'

generate_DepMap_models_ftINIT(chunk, use1Plus1)
