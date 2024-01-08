function [sol_refine, obj] = refine_solution(log)
%%
dof = double(log.dof);
Ns  = double(log.Ns);

t0 = log.param.theta_0;
xc_0 = log.param.xc_0;
yc_0 = 0.0;

lx = log.param.lx;
ly = log.param.ly;
m = log.param.mass;
I = log.param.inertial;
dt = log.param.control_param.dt;
g = -9.81;
%%
yalmip('clear')
% R = sdpvar(dof+1, 3, 3);
% F = sdpvar(dof+1, 3, 3);
% p = sdpvar(dof+1, 3, 1);
% v = sdpvar(dof+1, 3, 1);
% tau = sdpvar(dof, 3, 1);
% f = sdpvar(dof, 1);
% p_load = sdpvar(dof+1, 3, 1);
% v_load = sdpvar(dof+1, 3, 1);

x = sdpvar(dof+1, Ns+1); x_ = x';
y = sdpvar(dof,   Ns+1); y_ = y';
c = sdpvar(dof,   Ns+1); c_ = c';
s = sdpvar(dof,   Ns+1); s_ = s';

vx = sdpvar(dof+1, Ns+1); vx_ = vx';
vy = sdpvar(dof,   Ns+1); vy_ = vy';
ca = sdpvar(dof,   Ns+1); ca_ = ca';
sa = sdpvar(dof,   Ns+1); sa_ = sa';

lam_x = sdpvar(dof, Ns); lam_x_ = lam_x';
lam_y = sdpvar(dof, Ns); lam_y_ = lam_y';
tau   = sdpvar(Ns, 1);   tau_   = tau';

var = [x(:);  y(:);  c(:);  s(:);...
      vx(:); vy(:); ca(:); sa(:);...
      lam_x(:); lam_y(:); tau(:)];

% var = [x_(:);  y_(:);  c_(:);  s_(:);...
%       vx_(:); vy_(:); ca_(:); sa_(:);...
%       lam_x_(:); lam_y_(:); tau_(:)];

ns = dof * (Ns+1);
nf = Ns;
nlam = dof * Ns;

%%
% load
eq = [];
ineq = [];
f = 0;

eq = [eq;
    x(end, 1)  - xc_0;
    vx(end, 1) - 0];
%     c(1, 1) - cos(t0);
%     s(1, 1) - sin(t0);
%     cos(t0) * lx - sin(t0) * ly + x(1, 1) - xc_0;
%     sin(t0) * lx + cos(t0) * ly + y(1, 1) - yc_0;
%     vx(1, 1) - 0;
%     vy(1, 1) - 0;
%     ca(1, 1) - 1.0;
%     sa(1, 1) - 0];

for k = 1:dof
    eq = [eq;
        c(k, 1) - cos(t0(k));
        s(k, 1) - sin(t0(k));
        vx(k, 1) - 0.0;
        vy(k, 1) - 0.0;
        ca(k, 1) - 1.0;
        sa(k, 1) - 0.0];
    if k == 1
        eq = [eq;
            (cos(t0(k)) * lx - sin(t0(k)) * ly + x(k, 1) - x(end, 1) );
            (sin(t0(k)) * lx + cos(t0(k)) * ly + y(k, 1) - 0 )];
    else
        eq = [eq;
            (cos(t0(k)) * lx - sin(t0(k)) * ly + x(k, 1) - x(k-1, 1));
            (sin(t0(k)) * lx + cos(t0(k)) * ly + y(k, 1) - y(k-1, 1))];
    end
end



for k = 1:Ns
    %% cart
    eq = [eq;
        x(end, k+1) - (x(end, k) + dt * vx(end, k));
        (m * vx(end, k+1) - m * vx(end, k) - (lam_x(end, k) + tau(k))* dt )];

    for j = 1:dof
        eq = [eq;
            %% pole
            c(j, k+1) - (c(j, k) * ca(j, k) - s(j, k) * sa(j, k));
            s(j, k+1) - (s(j, k) * ca(j, k) + c(j, k) * sa(j, k));
            x(j, k+1) - (dt * (c(j, k) * vx(j, k) - s(j, k) * vy(j, k)) + x(j, k));
            y(j, k+1) - (dt * (s(j, k) * vx(j, k) + c(j, k) * vy(j, k)) + y(j, k))];
            %% sin cos
        ineq = [ineq; 
                -(s(j,k+1)^2   +  c(j,k+1)^2)  + 1.000000;
                  s(j,k+1)^2   +  c(j,k+1)^2   - 1 - 1e-8;
                -(sa(j,k+1)^2  +  ca(j,k+1)^2) + 1.000000;
                  sa(j,k+1)^2  +  ca(j,k+1)^2  - 1 - 1e-8;
                  x(j, k+1) + 1;
                 -x(j, k+1) + 1];
                %  tau(k) + 2;
                % -tau(k) + 2];

        %% kinematic constraints
        if j == 1
            eq = [eq;
                c(j,k+1) * lx - s(j,k+1) * ly + x(j, k+1) - x(end, k+1);
                s(j,k+1) * lx + c(j,k+1) * ly + y(j, k+1) - 0];
        else
            eq = [eq;
                c(j,k+1) * lx - s(j,k+1) * ly + x(j, k+1) - x(j-1, k+1);
                s(j,k+1) * lx + c(j,k+1) * ly + y(j, k+1) - y(j-1, k+1)];
        end

        Jx1 = (-c(j, k+1) * ly - s(j, k+1) * lx);
        Jy1 = (-s(j, k+1) * ly + c(j, k+1) * lx);
        eq = [eq; ((sa(j, k+1) - sa(j, k)) * I - (Jx1 * lam_x(j, k) + Jy1 * lam_y(j, k)) * dt^2)];

        Fx1 =   c(j, k+1)*lam_x(j, k) + s(j, k+1)*lam_y(j, k);
        Fy1 = - s(j, k+1)*lam_x(j, k) + c(j, k+1)*lam_y(j, k);

        if j < dof
            Fx2 =   c(j, k+1)*lam_x(j+1, k) + s(j+1, k+1)*lam_y(j, k);
            Fy2 = - s(j, k+1)*lam_x(j+1, k) + c(j+1, k+1)*lam_y(j, k);
        else
            Fx2 = 0.;
            Fy2 = 0.;
        end

        eq = [eq;
            (m * vx(j, k+1) - m * ( ca(j, k) * vx(j, k) + sa(j, k) * vy(j, k) ) - (Fx1 + Fx2 + m * s(1, k+1) * g)* dt );
            (m * vy(j, k+1) - m * (-sa(j, k) * vx(j, k) + ca(j, k) * vy(j, k) ) - (Fy1 + Fy2 + m * c(1, k+1) * g)* dt )];

        f = f + ( x(j, k)^2 + (y(j, k) - j * abs(ly))^2 +  (c(j, k) - 1)^2 +  s(j, k)^2)  * log.param.run_cost(1);
        f = f + (vx(j, k)^2 + vy(j, k)^2 + (ca(j, k) - 1)^2 + sa(j, k)^2) * log.param.run_cost(1);
        % f = f + (vx(j, k)^2 + vx(j, k)^2 + vy(j, k)^2   + (ca(j, k) - 1)^2 + sa(j, k)^2)  * log.param.run_cost(1);
        f = f + log.param.run_cost(1) * tau(k)^2;
    end
    f = f + ( x(end, k)^2 + vx(end, k)^2) * log.param.run_cost(1);
end
%%
f = f + ( x(end, end)^2 + vx(end, end)^2) * log.param.run_cost(1);
for j = 1:dof
    f = f + ( x(j, end)^2 + (y(j, end) - j * abs(ly))^2 +  (c(j, end) - 1)^2 +  s(j, end)^2) * log.param.run_cost(1);
    f = f + (vx(j, end)^2 + vx(j, end)^2 + vy(j, end)^2 + (ca(j, end) - 1)^2 + sa(j, end)^2) * log.param.run_cost(1);
end

for k = 1:length(var)-Ns
    ineq = [ineq;
        100 - var(k);
        100 + var(k)];
end

for k = 1:Ns+1
    ineq = [ineq;
        0 + ca(k)];
end

% ineq = [ineq; 
%         vx(:) + 5;
%        -vx(:) + 5;
%         vy(:) + 5;
%        -vy(:) + 5];
 
% for k = 1:Ns
%     ineq = [ineq;
%         5 - tau(k);
%         5 + tau(k);
%         55 - lam_x(k);
%         55 + lam_x(k);
%         55 - lam_y(k);
%         55 + lam_y(k);
%         ca(1, k)];
% end

% f = f + 1e5 * sum((var - log.sol_approx).^2);

% ineq = [ineq; 
%          var - (log.sol_approx - 6);
%         -var + (log.sol_approx + 6)];
%%
assign(var, log.sol_approx + 1e-2 * randn(size(log.sol_approx)))
c_norm = sqrt(value(c).^2 + value(s).^2);
c = c ./ c_norm;
s = s ./ c_norm;

ca_norm = sqrt(value(ca).^2 + value(sa).^2);
ca = ca ./ ca_norm;
sa = sa ./ ca_norm;

%%
ops = sdpsettings('solver','fmincon');%,'fmincon.algorithm','sqp');
% ops = sdpsettings('solver','ipopt');
ops.print_level = 12;
ops.verbose = true;

% % ops.fmincon.MaxIterations = 1e7;
% ops.ipopt.max_iter = 1; % 1e4 * 10;
% ops.ipopt.max_cpu_time = 1e100;

ops.usex0 = 1;
optimize([eq == 0; ineq >= 0],f,ops);
%%
x_pole = value(x(1, :));
x_cart = value(x(2, :));
y_pole = value(y(1, :));
c_pole = value(c(1, :));
s_pole = value(s(1, :));

vx_pole = value(vx(1, :));
vx_cart = value(vx(2, :));
vy_pole = value(vy(1, :));
ca_pole = value(ca(1, :));
sa_pole = value(sa(1, :));

lam_x = value(lam_x);
lam_y = value(lam_y);
tau   = value(tau);

obj = value(f);
%%
% sol_refine = [x_pole(:); x_cart(:); y_pole(:); c_pole(:); s_pole(:);...
%     vx_pole(:); vx_cart(:); vy_pole(:); ca_pole(:); sa_pole(:);...
%     lam_x(:); lam_y(:); tau(:)];
sol_refine = value(var);
end