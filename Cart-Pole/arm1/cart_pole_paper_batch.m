clear all;close all;clc;
%%
% dec_th3_3:   deg45
% th3_23:      deg90
% dec_th3_8:   deg135
% dec_th3_13:  deg171
%%
file_list = ["log45.mat",
             "log90.mat",
             "log135.mat",
             "log171.mat"];

% figure()
fig_size = [18 * 0.5, 1.3] * 6 ;
h = figure('Renderer', 'painters',  'unit', 'centimeters', 'Position', [0, 0, fig_size]);
% left_coner_list = [0.05, 0.3, 0.55, 0.8];

subplot(1, 1, 1, 'Position', [0.03, 0.01, 0.95, 1])

loc_offset = 0;
for k = 1:4

    load(file_list(k))
    
    if k == 1 || k == 2
        sol_eval= sol_approx;
        Ns = double(Ns);
        % sol_eval= sol;
        ns1 = Ns + 1;
        x_pole = sol_eval(1:2:ns1*2);
        x_cart = sol_eval(2:2:ns1*2);
        y_pole = sol_eval(ns1*2+1:3*ns1);
        c_pole = sol_eval(ns1*3+1:4*ns1);
        s_pole = sol_eval(ns1*4+1:5*ns1);

        vx_pole = sol_eval(5*ns1+1:2:ns1*7);
        vx_cart = sol_eval(5*ns1+2:2:ns1*7);
        vy_pole = sol_eval(7*ns1+1:ns1*8);
        ca_pole = sol_eval(ns1*8+1:9*ns1);
        sa_pole = sol_eval(ns1*9+1:10*ns1);

        offset = 10*ns1;
        lam_x = sol_eval(offset+1:offset+Ns);
        lam_y = sol_eval(offset+Ns+1:offset+2*Ns);
        tau   = sol_eval(offset+Ns*2 + 1:offset+3*Ns);
    else

    end

    x_pole = x_pole + loc_offset;
    x_cart = x_cart + loc_offset;

    
    hold on
    plot(x_cart, x_cart * 0, "bo")

    for j = 1:length(x_pole)
        haha = plot([x_pole(j), x_cart(j)], [y_pole(j), 0], "b-", "LineWidth", 1);
    end
    plot(x_pole, y_pole, "r-o", "LineWidth", 2)
    plot(x_pole(1), y_pole(1), "go", "LineWidth", 1, 'MarkerFaceColor', 'g')
    plot(x_pole(end), y_pole(end), "go", "LineWidth", 1, 'MarkerFaceColor', 'g')
    daspect([1, 1, 1])
    ylim([-0.5, 0.5])
    xlim([-0.5, 7.5])
    loc_offset = loc_offset + 2;

    % (obj - opt) / obj
end
grid on
box on
set(gca,'xtick',[])
set(h, 'Units','pixels');
set(h, 'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize', fig_size)
print("cp_all", "-dpdf")
%%
deglist = [45, 90, 135, 171];
radlist = [0.25, 0.5, 0.75, 0.95];

fig_size = [16, 12];
h = figure('Renderer', 'painters',  'unit', 'centimeters', 'Position', [0, 0, fig_size]);
for k = 1:4
    % hold on

    load(file_list(k))

    if k == 1 || k == 2
        sol_eval= sol_approx;
        Ns = double(Ns);
        % sol_eval= sol;
        ns1 = Ns + 1;
        x_pole = sol_eval(1:2:ns1*2);
        x_cart = sol_eval(2:2:ns1*2);
        y_pole = sol_eval(ns1*2+1:3*ns1);
        c_pole = sol_eval(ns1*3+1:4*ns1);
        s_pole = sol_eval(ns1*4+1:5*ns1);

        vx_pole = sol_eval(5*ns1+1:2:ns1*7);
        vx_cart = sol_eval(5*ns1+2:2:ns1*7);
        vy_pole = sol_eval(7*ns1+1:ns1*8);
        ca_pole = sol_eval(ns1*8+1:9*ns1);
        sa_pole = sol_eval(ns1*9+1:10*ns1);

        offset = 10*ns1;
        lam_x = sol_eval(offset+1:offset+Ns);
        lam_y = sol_eval(offset+Ns+1:offset+2*Ns);
        tau   = sol_eval(offset+Ns*2 + 1:offset+3*Ns);
    else

    end

    % x_pole = x_pole + loc_offset;
    % x_cart = x_cart + loc_offset;
    subplot(4, 1, k)
    % subplot(4, 1, k, 'Position', [0.05, (4 - k) * 0.26, 0.8, 0.2])
    k
    hold on
    plot(linspace(0, 5, 80), lam_x, "-o", "MarkerSize", 2)
    plot(linspace(0, 5, 80), lam_y, "-o", "MarkerSize", 2)
    plot(linspace(0, 5, 80), tau, "-o", "MarkerSize", 2)

    % title("$\theta_0 = " + num2str(deglist(k)) + " \deg$", "Interpreter","Latex")
    title("$\theta_0 = " + num2str(radlist(k)) + " \pi$", "Interpreter","Latex")

    % plot(lam_x, "b.", "MarkerSize", 5)
    % plot(lam_y, "b.", "MarkerSize", 5)
    % plot(tau, "b.", "MarkerSize", 5)
    % (obj - opt) / obj
    box on
    grid on
end
xlabel("Time (s)", "Interpreter","Latex")
legend({"$\lambda_x$", "$\lambda_y$", "$u$"}, "Interpreter","Latex")
%%
set(h, 'Units','pixels');
set(h, 'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize', fig_size)
print("cp_input", "-dpdf")