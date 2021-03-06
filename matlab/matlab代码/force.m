%% newmark方法的外载荷列向量，%每个向量元素是时间t的函数
%输入：t时间
%输入：NFreeDof是施加位移边界条件后，系统总的自由度的数目
%输出：f是系统的外载荷向量，即力边界条件
function f = force(t,NFreeDof)
f = zeros(NFreeDof,1);
%% 下面为每个自由度对应的外载荷赋予新值，默认的话就是零
%f(3) = 6;                          %书上例子的外载荷
%f(99) = 100*sin(2*pi*5*t);         %测试网上悬臂梁中点载荷的外载荷
%f(199)=0.8*sin(26.39*t);           %测试两节点的模型用
%f(596)=1;                          %测试自己悬臂梁梁的模型，100个单元，在末端加恒集中力外载荷

%f(596)=sin(31.4*t);                %测试自己悬臂梁梁的模型，100个单元，在末端加正弦载荷

% if(t==0)
%     f(596)=1;                     %测试自己悬臂梁梁的模型，100个单元，在末端加冲击
f(2)=sin(31.4*t);                %测试自己悬臂梁梁的模型，1个单元，在末端加正弦载荷
end