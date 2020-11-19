function createP90Figure(yvector1)
%CREATEFIGURE(yvector1)
%  YVECTOR1:  bar yvector

%  �� MATLAB �� 13-Apr-2019 12:47:04 �Զ�����

% ���� figure
figure1 = figure;

% ���� axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

% ���� bar
bar(yvector1,...
    'FaceColor',[0.635294117647059 0.07843137254902 0.184313725490196]);

% ���� ylabel
ylabel('90th percentile error','FontSize',20);

box(axes1,'on');
% ������������������
set(axes1,'XTick',[1 2],'XTickLabel',{'Scenario 1','Scenario 2'});
