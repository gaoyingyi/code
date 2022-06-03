clc;clear;
P=[0,0,0];
C=[1,1,1,1,0,0,0,0]';
Cf=[-eye(4);1,-1,0,0;0,-1,0,0;0,0,1,-1;-1,0,0,1];
q=[2 2 2 2 1 2 3 4]';
Q=diag(q);
fixedNodeCoordinates=[100,0,0;0,100,0;-100,0,0;0,-100,0;];
[x,y,z]=ForceDensityMethod(Q,C,Cf,P,fixedNodeCoordinates);