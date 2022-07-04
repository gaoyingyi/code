function F=beamInternalForce(deltaL)
E=350000;
r=5;
A=pi*r^2;
L=400;
F=E*A*deltaL/L;
end