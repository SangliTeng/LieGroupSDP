clear all;clc;close all;
%%
log = load("deg45.mat");
load("deg45_refine.mat")
% [sol_refine, obj] = refine_solution(log);
% save("deg45_refine.mat", "sol_refine")
%%
load("deg45.mat")
% figure()
sol_eval= sol_refine;
% sol_eval= sol;
ns1 = Ns + 1;

%%
fig_size = [7, 3] * 3 ;
h = figure('Renderer', 'painters',  'unit', 'centimeters', 'Position', [0, 0, fig_size]);
left_coner_list = [0.08, 0.55];

for id = 1:2
    subplot(1, 2, id, 'Position', [left_coner_list(id), 0.05, 0.4, 1]);
    hold on
    if id == 1
        [x, y, c, s, vx, vy, ca, sa, lam_x, lam_y, tau] = get_sol_eval(sol, ns1, dof, Ns);
        title("SDP solution", "Interpreter", "latex")
    else
        [x, y, c, s, vx, vy, ca, sa, lam_x, lam_y, tau] = get_sol_eval(sol_refine, ns1, dof, Ns);
        title("IPOPT refined solution", "Interpreter", "latex")
    end
    

    x_batch = reshape(x, [dof+1., ns1 + 0.]);
    y_batch = reshape(y, [dof+.0, ns1 + 0.]);
    for k = 1:length(x_batch)
        % plot([x_batch(3, k), x_batch(1, k)], [0, y_batch(1, k)], "b-", "LineWidth", 2)
        plot([x_batch(3, k), x_batch(1, k)], [0, y_batch(1, k)],             "-", "LineWidth", 1, "color", [1, 1, 1] * 0.25)
        plot([x_batch(2, k), x_batch(1, k)], [y_batch(2, k), y_batch(1, k)], "-", "LineWidth", 1, "color", [1, 1, 1] * 0.75)
    end

    plot(x(dof+1:dof+1:end), x(dof+1:dof+1:end) * 0, 'o', "LineWidth", 1, "color", [1, 1, 1] * 0.5, "MarkerSize", 3)
    for k = 2
        plot(x(k:dof+1:end), y(k:dof:end), 'r-o', "LineWidth", 1, "MarkerSize", 4)
        plot(x(k), y(k), 'go', "LineWidth", 1, "MarkerFaceColor", 'g')
        plot(x(end-1), y(end), 'go', "LineWidth", 1, "MarkerFaceColor", 'g')
    end

    xlabel("$x$", "Interpreter","latex")
    ylabel("$y$", "Interpreter","latex")

    grid on
    box on
    daspect([1, 1, 1])

    xlim([-1.2, 1.2])
    ylim([-0.8, 1.2])
end

set(h, 'Units','pixels');
set(h, 'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize', fig_size)
print("traj_cartpole3", "-dpdf")
%%
function [x, y, c, s, vx, vy, ca, sa, lam_x, lam_y, tau] = get_sol_eval(sol_eval, ns1, dof, Ns)
x = sol_eval(1:ns1*(dof+1));
y = sol_eval(ns1*(dof+1)+1:(dof * 2+1)*ns1);
c = sol_eval(ns1*(dof*2+1)+1:(dof*3+1)*ns1);
s = sol_eval(ns1*(dof*3+1)+1:(dof*4+1)*ns1);

offset = length([x; y; c; s]);
vx = sol_eval(offset + 1:offset + ns1*(dof+1));
vy = sol_eval(offset + ns1*(dof+1)+1:offset + (dof * 2+1)*ns1);
ca = sol_eval(offset + ns1*(dof*2+1)+1:offset + (dof*3+1)*ns1);
sa = sol_eval(offset + ns1*(dof*3+1)+1:offset + (dof*4+1)*ns1);
offset = length([x; y; c; s; ...
    vx; vy; ca; sa]);
% ns = dof * (Ns+1);
lam_x = sol_eval(offset + 1: offset + dof * Ns);
lam_y = sol_eval(offset + dof * Ns + 1: offset + dof * Ns * 2);
tau   = sol_eval(offset + dof * Ns * 2: end);
end
