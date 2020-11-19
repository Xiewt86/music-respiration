function [] = anim_plot(C, down, up)
figure
for i = 1 : length(C)
    plot(C{i})
    ylim([down, up])
    drawnow
end
end
