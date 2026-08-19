clear;
clc;

%% Load saved plot data

load('Figure10a_boxplot.mat');

selected_labels = plot_data.reaction_labels;

%% Prepare data for boxplot

plot_values = [];
box_groups = [];

for j = 1:numel(selected_labels)

    net_vals = plot_data.net{j};
    xch_vals = plot_data.xch{j};

    net_position = 2*j - 1;
    xch_position = 2*j;

    plot_values = [plot_values;
                   net_vals;
                   xch_vals];

    box_groups = [box_groups;
                  repmat(net_position, numel(net_vals), 1);
                  repmat(xch_position, numel(xch_vals), 1)];
end


%% Plot

figure;

boxplot(plot_values, box_groups, ...
    'Symbol','', ...
    'Positions',1:10);

hold on;


%% Colors

net_color = [0.2 0.45 0.75];
xch_color = [0.85 0.35 0.25];

boxes = findobj(gca,'Tag','Box');

for k = 1:numel(boxes)

    xdata = get(boxes(k),'XData');
    ydata = get(boxes(k),'YData');

    xpos = round(mean(xdata));

    if mod(xpos,2) == 1
        patch(xdata, ydata, net_color, ...
            'FaceAlpha',0.5, ...
            'EdgeColor',net_color);
    else
        patch(xdata, ydata, xch_color, ...
            'FaceAlpha',0.5, ...
            'EdgeColor',xch_color);
    end
end


%% Formatting

ylim([-1 2]);
xlim([0.5 10.5]);

xticks([1.5 3.5 5.5 7.5 9.5]);
xticklabels(selected_labels);

xlabel('Reaction');
ylabel('Flux');

title('Predicted Net and Exchange Fluxes');

box off;


%% Legend

hNet = patch(NaN,NaN,net_color, ...
    'FaceAlpha',0.5, ...
    'EdgeColor',net_color);

hXch = patch(NaN,NaN,xch_color, ...
    'FaceAlpha',0.5, ...
    'EdgeColor',xch_color);

legend([hNet hXch], {'Net','Exchange'}, ...
    'Location','best');