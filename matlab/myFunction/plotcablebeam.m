%% 根据单元节点的连接txt文件和节点坐标的txt文件绘制三维图
% elementFilename，单元文件：'ele.txt'
% NodeCoordinatesFilename,节点坐标文件：'nodeord.txt'
function []=plotcablebeam(elementFilename,nodeCoordinatesFilename)
    x=zeros(2,1);y=zeros(2,1);z=zeros(2,1);
    element=load(elementFilename);%加载文件
    [numElement,~]=size(element);%计算单元数量
    NodeCoordinates=load(nodeCoordinatesFilename);%加载文件
    [numNode,~]=size(NodeCoordinates);%计算节点数量
    hold on
    for j=1:numNode    %用于标记节点
        xnode=NodeCoordinates(j,2);
        ynode=NodeCoordinates(j,3);
        znode=NodeCoordinates(j,4);
        plot3(xnode,ynode,znode, ...
        'o','Color','blue','MarkerSize',10,'MarkerFaceColor','#FFFFFF')%修改节点的大小和颜色
        text(xnode,ynode,znode,num2str(j),'FontSize',16,'Color','blue')%修改节点标记的大小和颜色
    end
    
    for i=1:numElement  %用于绘制线段
        x(1,1)=NodeCoordinates(element(i,2),2);
        x(2,1)=NodeCoordinates(element(i,3),2);
        y(1,1)=NodeCoordinates(element(i,2),3);
        y(2,1)=NodeCoordinates(element(i,3),3);
        z(1,1)=NodeCoordinates(element(i,2),4);
        z(2,1)=NodeCoordinates(element(i,3),4);
        line(x,y,z,'LineWidth',3,'Color','black');%修改线宽和颜色
        text((x(1,1)+x(2,1))/2,(y(1,1)+y(2,1))/2,(z(1,1)+z(2,1))/2,num2str(i), ...
            'FontSize',20,'Color','red')%修改单元标记的大小和颜色
        view(3);
    end
end