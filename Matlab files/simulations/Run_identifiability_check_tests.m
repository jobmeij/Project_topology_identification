addpath(genpath('../toolbox')) % access the path with the functions



%% Example 3.11 from thesis
clear all;close all;clc

% This is the most straightforward example
net1.G = [0 1 1;1 0 1;1 1 0];
net1.R = [1 0;0 1;0 0]; 
net1.H = [1 1 0;1 1 0;0 0 1]; 


[identifiable,identifiable_nodes,identifiable_modules] = test_identifiability(net1) % test
% By manual check the network is identifiable




%% Example 3.21 from CDC2018 paper
clear all;close all;clc

% This example has no parameters in the noise model, and the noise is of a
% reduced rank, i.e. non-square H
net1.G = [0 0 0 0;0 0 0 0;0 1 0 0;1 1 1 0];
net1.R = []
net1.H = [eye(2);zeros(2)] ; 

fixed.H = net1.H; % The noise model has no parameters

[identifiable,identifiable_nodes,identifiable_modules] = test_identifiability(net1,fixed) % test
% By manual check node 4 is not identifiable





%% Example 3.14 from thesis (full G)
clear all;close all;clc

% This example has no parameters in R
net1.G = ones(5) - eye(5);
net1.R = [zeros(3,2);eye(2)];
net1.H = [1 1 0;1 1 0;0 0 1;zeros(2,3)] ;

fixed.R = net1.R; % has no parameters


[identifiable,identifiable_nodes,identifiable_modules] = test_identifiability(net1,fixed) % test
% By manual check the network is not identifiable, because of noise
% correlation in nodes 1 and 2.
% If we make noise 1 and 2 uncorrelated, then
net1.H = [1 0 0;0 1 0;0 0 1;zeros(2,3)] ;

[identifiable,identifiable_nodes,identifiable_modules] = test_identifiability(net1) % test


%% Example 3.17 from thesis (topology included)
clear all;close all;clc

% In this example the topology of G is included
net1.G = [0 1 0 0 1;1 0 1 0 0;0 0 0 1 0;0 1 0 0 0;0 0 1 0 0]; 
net1.R = [zeros(3,2);eye(2)];
net1.H = [1 1 0;1 1 0;0 0 1;zeros(2,3)] ; 

fixed.G = zeros(5); fixed.G(4,2) = 1; % G_42 is fixed and contains no parameters
fixed.R = net1.R; % has no parameters

[identifiable,identifiable_nodes,identifiable_modules] = test_identifiability(net1,fixed) % test
% By manual check the network is  identifiable




%% Example 3.4 from thesis
clear all;close all;clc

net1.G = [0 1 1;1 0 1;1 1 0];
net1.R = [1 0;0 1;1 0]; % The R matrix is fixed values
net1.H = []; % There is no noise model


fixed.R = net1.R; % The R matrix is fixed values


[identifiable,identifiable_nodes,identifiable_modules] = test_identifiability(net1,fixed) % test
% By manual check the network is indeed identifiable, in the generic sense
% However in the thesis we conclude that there are networks that are not
% identifiable. This is possible since we check generic network
% identifiability, and those non-identifiable networks form a set of
% measure 0.

