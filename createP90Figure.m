function createP90Figure(yvector1)
%CREATEFIGURE(yvector1)
%  YVECTOR1:  bar yvector

%  由 MATLAB 于 13-Apr-2019 12:47:04 自动生成

% 创建 figure
figure1 = figure;

% 创建 axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

% 创建 bar
bar(yvector1,...
    'FaceColor',[0.635294117647059 0.07843137254902 0.184313725490196]);

% 创建 ylabel
ylabel('90th percentile error','FontSize',20);

box(axes1,'on');
% 设置其余坐标区属性
set(axes1,'XTick',[1 2],'XTickLabel',{'Scenario 1','Scenario 2'});
