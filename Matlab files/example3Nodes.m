%% 3 node example for rank-reduced noise experiment
clear all; close all; clc;
% - 3 nodes w1, w2 and w3
% -- node w1 is connected to w2 through G12
% -- node w1 is connected to w3 through G13
% - 2 noise sources e1 and e2:
% -- e1 is connected to w1 through H11
% -- e2 is connected to w2 through H22 and to w3 through H23

% Create symbols
syms w1 w2 w3                   % nodes
syms G12 G13 G21 G23 G31 G32    % node transfer functions
syms H11 H22 H23                % noise transfer functions
syms e1 e2                      % noise sources
syms Qa lambda                  % for determining parameterized weighting matrix

% Create G matrix
G0 = [0 0 0; G12 0 0; G13 0 0]
G0_full = [0, G21, G31; G12, 0, G32; G13, G23, 0]

% Gamma = lim(z->inf){Hb0}, so if H23 (since Hb0 = [0 H23]) is monic, 
% Gamma = [0 1], else Gamma = [0 0]
Gamma = [0 1]   % It is stated that Hb0 must be monic, thus Gamma = [0 1]
% Gamma = [0 1 1; 0 0 1]

% Create H matrices
Ha0 = [H11, 0; 0, H22];
Hb0 = [0 H23];
% Hb0 = [0 H23 H23];
Hb0_ = Hb0 - Gamma
H0 = [Ha0; Hb0]
Hs0 = [H11 0 0; 0 H22 0; 0 (H23-1) 1]
Hd0 = [H11, 0 , 0; 0, H22, 0; 0, 0, eye(size(Hb0_,1))+Hb0_*pinv(Gamma)]

% nodes vector
w = [w1; w2; w3]
e = [e1; e2]
e_ = [e1; e2; e2]

% Identification criterion: w^= w-inv(Hs0)*((eye(3)-G)*w)
X_s = inv(Hs0)*(eye(3)-G0)
X_s_full = inv(Hs0)*(eye(3)-G0_full)
w_s_roof = w-X_s*w
w_s_roof_full = w-X_s_full*w
% all Hd0 related matrices
X_d = inv(Hd0)*(eye(3)-G0)
X_d_full = inv(Hd0)*(eye(3)-G0_full)
w_d_roof = w-X_d*w
w_d_roof_full = w-X_d_full*w
%
Q_theta = [Qa+lambda*Gamma'*Gamma, -lambda*Gamma'; -lambda*Gamma, lambda]
%
% epsilon = w-w_roof;
% epsilon_a = epsilon(1:2)
% epsilon_b = epsilon(3)
% Qaa = eye(2)

% Z_theta = (Gamma*epsilon_a-epsilon_b)'*(Gamma*epsilon_a-epsilon_b)

% argmin = epsilon_a'*Qaa*epsilon_a+lambda*Z_theta

% Prediction error
% epsilon_a = 
% PE = 
