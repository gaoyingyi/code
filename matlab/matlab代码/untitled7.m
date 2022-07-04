%% Write ABAQUS generated Matrix into Matlab
clc;clear
n=4;
ABAQUS_Matrix_rawdata=load('AAA_STIF2.mtx');

ABAQUS_Matrix=zeros(3*n,3*n);

for i=1:1:max(size(ABAQUS_Matrix_rawdata))

ABAQUS_Matrix(ABAQUS_Matrix_rawdata(i,1),ABAQUS_Matrix_rawdata(i,2))=ABAQUS_Matrix_rawdata(i,3);

end

%%%%%%%%%%%%%%%