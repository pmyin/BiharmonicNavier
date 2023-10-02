close all
clear
clc

format long

node=[0 0;0.5 0;0 0.5;0.5 0.5;1 0.5;0 1;0.5 1;1 1]-0.5;
node=node*4;
elem=[4 5 8; 4 8 7; 4 7 6; 4 6 3; 4 3 1; 4 1 2];
bdFlag = setboundary(node,elem,'Dirichlet');

mesh = struct('node',node,'elem',elem,'bdFlag',bdFlag); 

option.L0 = 0;
option.maxIt = 7;
option.printlevel = 1;
option.plotflag = 1;
option.elemType = 'P1';
option.tol = 1e-16;
option.fquadorder = 9;
option.quadorder = option.fquadorder;
pde = sincosdata;
femPoisson(mesh,pde,option);
