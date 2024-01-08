clear all;clc;close all;
%%
dof_list = [6, 7, 10, 15, 20, 25, 30, 40];
id = 50;
%%
ws_list = {};
elapsed_list = [];
for d = 1:length(dof_list)
    for k = 1:id
        file_name = "log_" + num2str(dof_list(d)) + "_" + num2str(k);
        ws = load(file_name);
        ws_list{d, k}  = ws;
        elapsed_list(d, k) = ws.elapsed;
    end
end
%%
boxplot(elapsed_list', "Positions", dof_list, "Labels", dof_list, "Colors", "b", "Widths", 1, 'Whisker',1000)
box on
grid on
hold on
plot(dof_list, median(elapsed_list'), "k--")

xlabel("DOF of Manipulator", "fontsize", 12)
ylabel("Time (s)", "fontsize", 12)