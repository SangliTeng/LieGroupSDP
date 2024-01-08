clear;clc;close all;
%%
case_num = 5;
load("log_NC_" + num2str(case_num) + ".mat")
%%
dof = double(dof);
dx0 = double(dx0);
dz0 = double(dz0);
%% warm start
% [R, F, p, v, tau, f, p_load, v_load, lam, opt_fmincon] = refine_solution_det_load(sol + randn(size(sol)) * 0e-2, dof, dt, Mass, Inertial, Qc, Rc, Pc, double(quat_init'), obs_flag, dx0, dz0);
% [R, F, p, v, tau, f, opt_fmincon] = refine_solution(sol + randn(size(sol)) * 0e-4, dof, dt, Mass, Inertial, Qc, Rc, Pc, double(quat_init'), obs_flag);
% load("ipopt-refine/refine_NC_23.mat")
% load("ipopt-refine3/refine_NC_" + num2str(case_num) + ".mat")
load("refine_NC_" + num2str(case_num) + ".mat")

%%
R_tssos = sol(1:9 * (dof + 1));
F_tssos = sol(1+9 * (dof + 1):9 * (dof + 1) * 2);

shift = 9 * (dof + 1) * 2;
p_tssos = sol(          1+shift : shift+3*(dof+1));
v_tssos = sol(1+shift+3*(dof+1) : shift+6*(dof+1));
v_tssos = reshape(v_tssos, [dof+1, 3]);
p_tssos = reshape(p_tssos, [dof+1, 3]);

shift = shift+6*(dof+1);
control = sol(shift+1:shift + 4 * dof);

shift = shift+4*dof;
p_tssos_load = sol(          1+shift : shift+3*(dof+1));
v_tssos_load = sol(1+shift+3*(dof+1) : shift+6*(dof+1));
v_tssos_load = reshape(v_tssos_load, [dof+1, 3]);
p_tssos_load = reshape(p_tssos_load, [dof+1, 3]);


shift = shift + 6 * (dof + 1);
lam_tssos = sol(1+shift : end);

R_tssos = reshape(R_tssos, [dof+1, 3, 3]);
F_tssos = reshape(F_tssos, [dof+1, 3, 3]);


tau_tssos = reshape(control(1:dof*3), [dof, 3]);
f_tssos = reshape(control(dof*3+1:end), [dof, 1]);
%%
% eq_tssos = test_feasibility(R_tssos, F_tssos, p_tssos, v_tssos, tau_tssos, f_tssos, dof);
% eq_fmincon = test_feasibility(R, F, p, v, tau, f, dof);

%%
w = zeros(dof+1, 1);
for k = 1:dof + 1
    FF = squeeze(F(k,:,:));
    w(k) = sqrt(sum(logm(FF).^2, "all") / 2);
end
%%
% subplot(1,2,1)
x = p(1:dof+1);
x_load = p_load(1:dof+1);
% x_load = p_tssos_load(1:dof+1);
% dx = ones(dof+1, 1) * 0.1;
y = p(dof+1+1:2*(dof+1));
y_load = p_load(dof+1+1:2*(dof+1));
% y_load = p_tssos_load(dof+1+1:2*(dof+1));
% dy = zeros(dof+1, 1) * 0.1;
z = p(2*(dof+1)+1:end);
z_load = p_load(2*(dof+1)+1:end);
% z_load = p_tssos_load(2*(dof+1)+1:end);
% dz = zeros(dof+1, 1) * 0.1;

% figure()
fig_size = [18 * 0.5, 4.0] * 3 ;
h = figure('Renderer', 'painters',  'unit', 'centimeters', 'Position', [0, 0, fig_size]);

left_coner_list = [0.05, 0.38, 0.68];
% left_coner_list = [0.1, 0.6];

for fig_id = 1:2
    % subplot(1, 2, fig_id)
    subplot(1, 3, fig_id, 'Position', [left_coner_list(fig_id), 0.01, 0.3, 1]);
    hold on
    len = 0.1;
    for k = 1:dof + 1
        Rk = squeeze(R(k, :, :));
        dx = Rk(:, 1) * len;

        if mod(k, 3) == 1 || k <= 10 % || true
            plot3([x(k), x(k) + dx(1)], [y(k), y(k) + dx(2)], [z(k), z(k) + dx(3)], "r-", "LineWidth",1);

            dy = Rk(:, 2) * len;
            plot3([x(k), x(k) + dy(1)], [y(k), y(k) + dy(2)], [z(k), z(k) + dy(3)], "g-", "LineWidth",1);

            dz = Rk(:, 3) * len;
            plot3([x(k), x(k) + dz(1)], [y(k), y(k) + dz(2)], [z(k), z(k) + dz(3)], "b-", "LineWidth",1);

            plot3([x(k), x_load(k)], [y(k), y_load(k)], [z(k), z_load(k)], "b--");
        end
    end
    plot3(x, y, z, "-")
    plot3(x_load, y_load, z_load, ".", "MarkerSize", 6)
    plot3(x_load, y_load, z_load)

    [X, Y] = meshgrid(-2:0.1:2, -2:0.1:2);
    fig1 = mesh(X, Y, X * 0);
    [X, Y, Z] = cylinder(0.5, 200);
    surf(X + 0.6, Y + 0.5, Z * 6, 'facecolor','r','LineStyle','none','facealpha',0.2)

    grid on
    box on
    xlim([-0.2 - 0.5, 1.5 + 0.5])
    ylim([-0.2 - 0.5, 1.5 + 0.5])
    zlim([-0.1, 3.5])
    daspect([1,1,1])

    if fig_id == 1
        view(45, 35)
    end

    if fig_id == 2
        view(-60, 70)
    end
    
    lw_pos = 1.2;
    len_arrow = 0.1;
    font_size = 14;
    legend_size = 8;

    xlabel("$x$", "Interpreter", "latex", "FontSize", font_size)
    ylabel("$y$", "Interpreter", "latex", "FontSize", font_size)
    zlabel("$z$", "Interpreter", "latex", "FontSize", font_size)
end
%%
fig_id = 3;
subplot(1, 3, fig_id, 'Position', [0.75, 0.1, 0.20, 0.8]);
hold on
plot(x_load, y_load, "r-o", "MarkerSize", 2)
plot(x, y, "b-o", "MarkerSize", 2)

plot(p_tssos_load(1:dof+1), p_tssos_load(dof+2:2*(dof + 1)), "r--", "MarkerSize", 2)
plot(p_tssos(1:dof+1), p_tssos(dof+2:2*(dof + 1)), "b--", "MarkerSize", 2)

% plot(x, y, "b.")
% plot(x_load, y_load, "b.")
tt = 0:0.01:2 * pi;
plot(cos(tt) * 0.5 + 0.6, 0.5 * sin(tt) + 0.5, "k--", "LineWidth", 1)
grid on
box on
daspect([1, 1, 1])
xlim([-0.1, 1.2])
ylim([-0.1, 1.6])
xlabel("$x$", "Interpreter", "latex")
ylabel("$y$", "Interpreter", "latex")
legend("$p^{L}$ - Refined", "$p$ - Refined", "$p^L$ - SDP", "$p$ - SDP", "location", "northwest", "Interpreter", "latex")

%%
% fig_id = 3;
% subplot(1, 3, fig_id, 'Position', [0.75, 0.1, 0.20, 0.8]);
% hold on
% plot(x_load, y_load, "-o", "MarkerSize", 2)
% plot(x, y, "-o", "MarkerSize", 2)
% % plot(x, y, "b.")
% % plot(x_load, y_load, "b.")
% tt = 0:0.01:2 * pi;
% plot(cos(tt) * 0.5 + 0.6, 0.5 * sin(tt) + 0.5, "k--", "LineWidth", 1)
% grid on
% box on
% daspect([1, 1, 1])
% xlim([-0.1, 1.2])
% ylim([-0.1, 1.6])
% xlabel("$x$", "Interpreter", "latex")
% ylabel("$y$", "Interpreter", "latex")
% legend("Load position", "Drone position", "location", "northwest", "Interpreter", "latex")
%%
set(h, 'Units','pixels');
set(h, 'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize', fig_size)
print("v2_load_traj_" + num2str(case_num), "-dpdf")

%%
(opt_fmincon - opt) / opt_fmincon