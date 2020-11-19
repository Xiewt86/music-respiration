function createNrmseFigure(yvector1)
%CREATEFIGURE(yvector1)
%  YVECTOR1:  bar yvector

%  由 MATLAB 于 13-Apr-2019 12:48:30 自动生成

% 创建 figure
figure('OuterPosition',...
    [353 184.333333333333 574.666666666667 506.666666666667]);

% 创建 axes
axes1 = axes;
hold(axes1,'on');

% 创建 bar
bar(yvector1,...
    'FaceColor',[0.635294117647059 0.07843137254902 0.184313725490196]);

% 创建 ylabel
ylabel('NRMSE','FontSize',20);

box(axes1,'on');
% 设置其余坐标区属性
set(axes1,'XTick',[1 2],'XTickLabel',{'Scenario 1','Scenario 2'});
