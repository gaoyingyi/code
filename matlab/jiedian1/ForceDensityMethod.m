%最原始的力密度法
% 通过力密度法输出自由节点的x,y,z三个方向的坐标
%有n个节点，m个单元,其中自由节点nd个，固定节点nf个，n=nd+nf
%已经经过简单测试，计算结果准确。
%     Q=[];%对角力密度矩阵，矩阵维度m×m，容易给出，等力密度准则都是常值
%     C=[];%自由节点拓扑矩阵，矩阵维度m×nd，需要代码生成
%     Cf=[];%固定节点拓扑矩阵，矩阵维度m×nf 需要代码生成
%     P=[];%自由节点外载荷，矩阵维度nd×3  容易给出，正常的话都是0
%     fixedNodeCoordinates=[];%固定节点坐标，矩阵维度nf×3，容易给出
function [x,y,z]=ForceDensityMethod(Q,C,Cf,P,fixedNodeCoordinates)
    x=(P(:,1)-C'*Q*Cf*fixedNodeCoordinates(:,1))/(C'*Q*C);
    y=(P(:,2)-C'*Q*Cf*fixedNodeCoordinates(:,2))/(C'*Q*C);
    z=(P(:,3)-C'*Q*Cf*fixedNodeCoordinates(:,3))/(C'*Q*C);
end