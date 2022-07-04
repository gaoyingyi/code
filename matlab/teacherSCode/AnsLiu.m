clc;
clear;
load Coordinate.txt;
load Beam.txt;
load Beam_XY.txt;
load Beam_Z.txt
load Force.txt
X0=1e-3*Coordinate(:,2);
Y0=1e-3*Coordinate(:,3);
Z0=1e-3*Coordinate(:,4);
[N,N1]=size(Coordinate);
[M_B,M_B1]=size(Beam);
[M_XY,M_XY1]=size(Beam_XY);
[M_Z,M_Z1]=size(Beam_Z);
[N_F,MF1]=size(Force);

fid=fopen('BeamJssAns.mac','wt');
fprintf(fid,'%s\n','FINISH');
fprintf(fid,'%s\n','/CLEAR');
fprintf(fid,'%s\n','*SET,ContactEx,3.5e11          !接触梁单元弹性模量');
fprintf(fid,'%s\n','*SET,ContactR,10e-6           !接触梁单元单元横截面积');
fprintf(fid,'%s\n','*SET,BeamEx,3.5e11         !梁单元弹性模量');
fprintf(fid,'%s\n','*SET,BeamPrxy,0.28         !梁单元泊松比');
fprintf(fid,'%s\n','*SET,BeamR,26e-6           !梁单元横截半径');
fprintf(fid,'%s\n','*SET,DF,0.001              !四周力');
fprintf(fid,'%s\n\n','*SET,PI,3.141592653');

fprintf(fid,'%s\n','/FILNAME,JSS');
fprintf(fid,'%s\n','/PREP7');
fprintf(fid,'%s\n','ET,1,BEAM188              !3-D 2-NODE BEAM');
fprintf(fid,'%s\n','KEYOPT,1,3,2');
fprintf(fid,'%s\n','MP,EX,1,BeamEx             !梁单元材料参数');
fprintf(fid,'%s\n','MP,PRXY,1,BeamPrxy ');
fprintf(fid,'%s\n','MP,DENS,1,19.35E3');
fprintf(fid,'%s\n','SECTYPE,1, BEAM, CSOLID');
fprintf(fid,'%s\n','SECDATA,BeamR');

fprintf(fid,'%s\n','ET,2,Beam188                !暂时不设置只受压');
fprintf(fid,'%s\n','KEYOPT,2,3,2');
fprintf(fid,'%s\n','MP,EX,2,ContactEx             !梁单元材料参数');
fprintf(fid,'%s\n','MP,PRXY,2,BeamPrxy ');
fprintf(fid,'%s\n','MP,DENS,2,19.35E3 ');
fprintf(fid,'%s\n','SECTYPE,2, BEAM, CSOLID');
fprintf(fid,'%s\n','SECDATA,ContactR');



for i=1:N
    fprintf(fid,'%s%f%s%f%s%f%s%f\n','N,',Coordinate(i,1),',',X0(i),',',Y0(i),',',Z0(i));
end
fprintf(fid,'\n%s\n','TYPE,1');
fprintf(fid,'%s\n','MAT,1');
fprintf(fid,'%s\n','SECNUM,1');
for i=1:M_B
    fprintf(fid,'%s%d%s%d%s%d\n','EN,',Beam(i,1),',',Beam(i,2),',',Beam(i,3));
end
fprintf(fid,'\n%s\n','TYPE,2');
fprintf(fid,'%s\n','MAT,2');
fprintf(fid,'%s\n','REAL,2');
fprintf(fid,'%s\n','SECNUM,2');
for i=1:M_XY
    fprintf(fid,'%s%d%s%d%s%d\n','EN,',Beam_XY(i,1),',',Beam_XY(i,2),',',Beam_XY(i,3));
end
fprintf(fid,'\n');
for i=1:M_Z
    fprintf(fid,'%s%d%s%d%s%d\n','EN,',Beam_Z(i,1),',',Beam_Z(i,2),',',Beam_Z(i,3));
end
for i=1:N_F
    fprintf(fid,'%s%d%s%f\n','F,',Force(i,1),',Fx,',Force(i,2));
    fprintf(fid,'%s%d%s%f\n','F,',Force(i,1),',Fy,',Force(i,3));
    fprintf(fid,'%s%d%s%f\n','F,',Force(i,1),',Fz,',Force(i,4));
end
%边界只留一个自由度
 for i=1:12
    fprintf(fid,'%s%d%s\n','D,',Force(i,1),',UX,,,,,UZ,ROTX,ROTY,ROTZ');
end
for i=13:20
    fprintf(fid,'%s%d%s\n','D,',Force(i,1),',UY,,,,,UZ,ROTX,ROTY,ROTZ');
end

%边界留两个自由度
% for i=1:N_F
%      fprintf(fid,'%s%d%s\n','D,',Force(i,1),',UZ,,,,,ROTX,ROTY,ROTZ');
% end

fprintf(fid,'\n%s\n','FINISH');

fprintf(fid,'\n%s\n','/SOLU');
fprintf(fid,'%s\n','ANTYPE,STATIC');
fprintf(fid,'%s\n','NLGEOM,ON');
fprintf(fid,'%s\n','TIMINT,ON');
fprintf(fid,'%s\n','OUTERS,ALL,ALL');
fprintf(fid,'%s\n','NSUBST,50,100,50');
fprintf(fid,'%s\n','NROPT,FULL');
fprintf(fid,'%s\n','KBC,1');
fprintf(fid,'%s\n','TIME,1');
fprintf(fid,'%s\n','ALLSEL,ALL');

fprintf(fid,'%s\n','SOLVE');
fprintf(fid,'%s\n','');

fclose('all');