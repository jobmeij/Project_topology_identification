%% new example
clear all; close all; clc

syms w1 w2 w3 w4 w5;
syms e1 e2 e3;
syms H11 H22 H33 H34 H45 H15 H42 H43 H51 H53;
syms G12 G15 G23 G54;

Hb = [0 H42 H43; H51 0 H53]

Gamma = [0 1 1; 1 0 1]

Hb_ = Hb-Gamma



w = [w1 w2 w3 w4 w5]';
e = [e1 e2 e3]';