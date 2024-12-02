clear;close all;clc;
load("logger_3.mat")
%%
sv_list = [];
for k = 1:length(moment)
    sv = svd(moment{k});
    sv_list = [sv_list;  sv(2) / sv(1)];
end

figure()
plot( -sort(-log10(sv_list)))

% sol = sol_approx;
%%
% Ns = dof;
sol_exec = sol;
% load("freespace.mat")
% [sol_refine, obj] = refine_solution_load_cvx_obs(sol, double(Ns), dt, 0.5, diag([1, 2, 1]) / 10, [1, 0, 0, 0]', A_list, b_list)
%%
figure()
subplot(2, 1, 1)
q = sol_exec(1:4*(Ns+1));
q = reshape(q, [Ns+1, 4]);
hold on
plot(q)

w = sol_exec(1+4*(Ns+1):2*4*(Ns+1));
w = reshape(w, [Ns+1, 4]) * dt / 2;
plot(w)

p = sol_exec(2*4*(Ns+1)+1:2*4*(Ns+1) + 3 * (Ns+1));
p = reshape(p, [Ns+1, 3]);

v = sol_exec(2*4*(Ns+1) + 3 * (Ns+1)+1:2*4*(Ns+1) + 2*3*(Ns+1));
v = reshape(v, [Ns+1, 3]);

p_load = sol_exec(2*4*(Ns+1) + 2*3*(Ns+1)+1:2*4*(Ns+1) + 3*3*(Ns+1));
p_load = reshape(p_load, [Ns+1, 3]);

%%
vertex_obs = {
    [[3.4, 2.6]; [3.4, 4.6];[2.4, 4.6]; [2.4, 2.6]; [1.4, 2.2]; [3.8, 0.2];[4.8, 1.2]],
    [[1.4, 2.8]; [2.2-0.2, 2.8]; [2.2-0.2, 4.6]; [1.4, 4.6]],
    [[1.0, 2.6+0.2]; [1.0, 5.0]; [0.4, 5.0]; [0.4, 2.6+0.2]],
    [[1.0, 2.4]; [1.0, 0.0]; [0.4, 0.0]; [0.4, 2.4]],
    [[3.8, 3.0+0.2]; [3.8, 5.0]; [4.4, 5.0]; [4.4, 3.0+0.2]],
    [[3.8, 2.8]; [3.8, 2.6]; [5.0, 2.6]; [5.0, 2.8]]};

vertex_free = {
    [[0.4, 0.0]; [0.4, 5.0]; [0.0, 5.0]; [0.0, 0.0]],
    [[0.4, 2.4]; [1.0, 2.4]; [1.0, 2.8]; [0.4, 2.8]],
    [[1.4, 2.2]; [1.4, 4.6]; [1.0, 4.6]; [1.0, 2.2]],
    [[1.4, 2.2]; [2.4, 2.6]; [2.4, 2.8]; [1.4, 2.8]],
    [[2.0, 2.8]; [2.4, 2.8]; [2.4, 4.6]; [2.0, 4.6]],
    [[1.4, 2.2]; [1.0, 2.2]; [1.0, 0.0]; [3.8, 0.0]; [3.8, 0.2]],
    [[3.8, 4.6]; [3.8, 5.0]; [1.0, 5.0]; [1.0, 4.6]],
    [[5.0, 0.0]; [5.0, 1.2]; [4.8, 1.2]; [3.8, 0.2]; [3.8, 0.0]],
    [[3.4, 2.6]; [4.8, 1.2]; [5.0, 1.2]; [5.0, 2.6]],
    [[3.4, 2.6]; [3.8, 2.6]; [3.8, 4.6]; [3.4, 4.6]],
    [[3.8, 2.8]; [4.4, 2.8]; [4.4, 3.2]; [3.8, 3.2]],
    [[5.0, 2.8]; [5.0, 5.0]; [4.4, 5.0]; [4.4, 2.8]],
    };
figure()
fig_size = [2, 0.6] * 20 ;
h = figure('Renderer', 'painters',  'unit', 'centimeters', 'Position', [5, 5, fig_size]);

for fig_id = 1:2

    if fig_id == 1
        offset = 0.2;
        line_style = '-';
        subplot(1, 2, 1, 'Position', [0.1, 0.1, 0.4, 0.9])
    else
        offset = 2;
        line_style = '-';
        subplot(1, 2, 2, 'Position', [0.1 + 0.5, 0.15, 0.3, 0.8])
    end
    % subplot(1, 2, 1, 'Position', [0.1, 0.15, 0.85, 0.8])

    % subplot(1, 2, fig_id)

    plot3(p(:, 1), p(:, 2), p(:, 3) + offset, '-o', "LineWidth", 1)
    hold on
    % for k = 1:length(q)
    %     % if mod(k, 3) == 0 || k <= 10 || k == logger{id}.dof + 1
    %         len_arrow = 0.1;
    %         lw_pos = 0.2;
    %         Rk = quat2rotm(q(k, :));
    %         dx = Rk(:, 1) * len_arrow;
    %         % plot3([x(k), x(k) + dx(1)], [p(k, 2), p(k, 2) + dx(2)], [p(k, 3), p(k, 3) + dx(3)], "r-", "LineWidth",lw_pos);
    %         plot3([p(k, 1) - dx(1),   p(k, 1) + dx(1)], [p(k, 2) - dx(2), p(k, 2) + dx(2)], [p(k, 3) - dx(3), p(k, 3) + dx(3)], "b-o", "LineWidth",lw_pos);
    %         plot3([p(k, 1) + dx(1)], [ p(k, 2) + dx(2)], [ p(k, 3) + dx(3)], "r-o", "LineWidth",lw_pos);
    %
    %         dy = Rk(:, 2) * len_arrow;
    %         % plot3([x(k), x(k) + dy(1)], [p(k, 2), p(k, 2) + dy(2)], [p(k, 3), p(k, 3) + dy(3)], "g-", "LineWidth",lw_pos);
    %         plot3([p(k, 1) - dy(1), p(k, 1) + dy(1)], [p(k, 2) - dy(2), p(k, 2) + dy(2)], [p(k, 3) - dy(3), p(k, 3) + dy(3)], "b-o", "LineWidth",lw_pos);
    %
    %         dz = Rk(:, 3) * len_arrow;
    %         % plot3([p(k, 1), p(k, 1) + dz(1)], [p(k, 2), p(k, 2) + dz(2)], [p(k, 3), p(k, 3) + dz(3)], "b-", "LineWidth",lw_pos);
    %         plot3([p(k, 1), p(k, 1) + dz(1)], [p(k, 2), p(k, 2) + dz(2)], [p(k, 3), p(k, 3) + dz(3)], "k-", "LineWidth",lw_pos);
    %     % else
    %
    %     % end
    % end

    hold on
    plot3(p_load(:, 1), p_load(:, 2), p_load(:, 3) + offset, '-o', "LineWidth", 1)
    daspect([1, 1, 1])

    for k = 1:length(p)-1
        xv = [p(k, 1), p_load(k, 1)];
        yv = [p(k, 2), p_load(k, 2)];
        zv = [p(k, 3), p_load(k, 3)] + offset;

        plot3(xv, yv, zv, 'k--')
    end
    % (p(:, 1) - p_load(:, 1)).^2 + (p(:, 2) - p_load(:, 2)).^2 + (p(:, 3) - p_load(:, 3)).^2;

    vertex_exec = vertex_obs;
    for k = 1:length(vertex_exec)
        v = vertex_exec{k};

        X_top = [];
        Y_top = [];
        for j = 1:size(v, 1)
            v1 = v(j, :);
            if j == size(v, 1)
                v2 = v(1, :);
            else
                v2 = v(j+1, :);
            end

            if fig_id == 2
                plot3([v1(1), v2(1)], [v1(2), v2(2)], [-1, -1], "k-")
                plot3([v1(1), v2(1)], [v1(2), v2(2)], [1, 1] * 0, "k-")
            end

            X = linspace(v1(1), v2(1), 2);
            X = [X(:)'; X(:)'];
            Y = linspace(v1(2), v2(2), 2);
            Y = [Y(:)'; Y(:)'];

            Z = [0, 0; 1, 1];
            surf(X, Y, Z - 1, 'facecolor','r','LineStyle',line_style,'facealpha',0.1)

            hold on

        end

        pv = polyshape(vertex_obs{k});
        plot(pv, 'facecolor','r','LineStyle',line_style,'facealpha',0.1)

        grid on
        box on
        xlim([0, 5])
        ylim([0, 5])

        xticks(0:1:5)
        yticks(0:1:5)
    
    end
    daspect([1, 1, 1])

    if fig_id == 1
        view(45, 45 - 10)
    elseif fig_id == 2
        view(0, 90)

        for k = [1, 2, 3, 4, 5, 7, 10, 11, 12] % 1:length(vertex_free)
            pv = polyshape(vertex_free{k});
            plot(pv, 'facecolor','g', 'LineStyle',line_style,'facealpha',0.1) %
        end

        for k = [6, 8, 9] % 1:length(vertex_free)
            pv = polyshape(vertex_free{k});
            plot(pv, 'facecolor',[1, 1, 1] / 2, 'LineStyle',line_style,'facealpha',0.1) %
        end
    end
end
set(h, 'Units','pixels');
set(h, 'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize', fig_size)
print("quat-obs", "-dpdf")
%%
figure()
lam = sol_exec(2*4*(Ns+1) + 4*3*(Ns+1)+1:2*4*(Ns+1) + 4*3*(Ns+1) + Ns);

%
tau = sol_exec(2*4*(Ns+1) + 4*3*(Ns+1) + Ns + 1:2*4*(Ns+1) + 4*3*(Ns+1) + Ns + 3*Ns);
tau = reshape(tau, [Ns, 3]) ;

plot(tau)
%
fk = sol_exec(2*4*(Ns+1) + 4*3*(Ns+1) + Ns + 3*Ns + 1:2*4*(Ns+1) + 4*3*(Ns+1) + Ns + 6*Ns);
fk = reshape(fk, [Ns, 3]);

plot(fk)
%%
r = sol_exec(end-Ns+1:end);
plot(r)