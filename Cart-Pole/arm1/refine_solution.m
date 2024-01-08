function [x_pole, x_cart, y_pole, c_pole, s_pole, vx_pole, vx_cart, vy_pole, ca_pole, sa_pole, lam_x, lam_y, tau, obj, sol_refine] = refine_solution(log)
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

x = sdpvar(dof+1, Ns+1);
y = sdpvar(dof,   Ns+1);
c = sdpvar(dof,   Ns+1);
s = sdpvar(dof,   Ns+1);

vx = sdpvar(dof+1, Ns+1);
vy = sdpvar(dof,   Ns+1);
ca = sdpvar(dof,   Ns+1);
sa = sdpvar(dof,   Ns+1);

lam_x = sdpvar(dof, Ns);
lam_y = sdpvar(dof, Ns);
tau   = sdpvar(Ns, 1);

ns = dof * (Ns+1);
nf = Ns;
nlam = dof * Ns;

%%
% load
eq = [];
ineq = [];
f = 0;

eq = [eq;
    x(2, 1)  - xc_0;
    vx(2, 1) - 0;
    c(1, 1) - cos(t0);
    s(1, 1) - sin(t0);
    cos(t0) * lx - sin(t0) * ly + x(1, 1) - xc_0;
    sin(t0) * lx + cos(t0) * ly + y(1, 1) - yc_0;
    vx(1, 1) - 0;
    vy(1, 1) - 0;
    ca(1, 1) - 1.0;
    sa(1, 1) - 0];


for k = 1:Ns
    Jx = (-c(1, k+1) * ly - s(1, k+1) * lx);
    Jy = (-s(1, k+1) * ly + c(1, k+1) * lx);
    eq = [eq;
        %% pole
        c(1, k+1) - (c(1, k) * ca(1, k) - s(1, k) * sa(1, k));
        s(1, k+1) - (s(1, k) * ca(1, k) + c(1, k) * sa(1, k));
        x(1, k+1) - (dt * (c(1, k) * vx(1, k) - s(1, k) * vy(1, k)) + x(1, k));
        y(1, k+1) - (dt * (s(1, k) * vx(1, k) + c(1, k) * vy(1, k)) + y(1, k));
        %% cart
        x(2, k+1) - (x(2, k) + dt * vx(2, k));
        %% kinematic constraints
        c(1,k+1) * lx - s(1,k+1) * ly + x(1, k+1) - x(2, k+1);
        s(1,k+1) * lx + c(1,k+1) * ly + y(1, k+1) - 0;
        s(1,k+1)^2  +  c(1,k+1)^2 - 1.0;
        sa(1,k+1)^2 + ca(1,k+1)^2 - 1.0;
        %% pole
        ((sa(1, k+1) - sa(1, k)) * I - (Jx * lam_x(1, k) + Jy * lam_y(1, k)) * dt^2);
        (m * vx(1, k+1) - m * ( ca(1, k) * vx(1, k) + sa(1, k) * vy(1, k) ) - (  c(1, k+1)*lam_x(1, k) + s(1, k+1)*lam_y(1, k) + m * s(1, k+1) * g)* dt );
        (m * vy(1, k+1) - m * (-sa(1, k) * vx(1, k) + ca(1, k) * vy(1, k) ) - (- s(1, k+1)*lam_x(1, k) + c(1, k+1)*lam_y(1, k) + m * c(1, k+1) * g)* dt );
        (m * vx(2, k+1) - m * vx(2, k) - (lam_x(1, k) + tau(k))* dt );
        ];

    f = f + ( x(2, k)^2 +  x(1, k)^2 + (y(1, k) - abs(ly))^2 + (c(1, k) - 1)^2 + s(1, k)^2 ) * log.param.run_cost(1);
    f = f + (vx(2, k)^2 + vx(1, k)^2 + vy(1, k)^2 + (ca(1, k) - 1)^2 + sa(1, k)^2) * log.param.run_cost(1);
    f = f + log.param.run_cost(1) * tau(k)^2;
end
%%
f = f + ( x(2, end)^2 +  x(1, end)^2 + (y(1, end) - abs(ly))^2 + (c(1, end) - 1)^2 + s(1, end)^2) * log.param.terminal_cost(1);
f = f + (vx(2, end)^2 + vx(1, end)^2 + vy(1, end)^2 + (ca(1, end) - 1)^2 + sa(1, end)^2) * log.param.terminal_cost(1);

for k = 1:Ns
    ineq = [ineq;
            20 - tau(k);
            20 + tau(k);
            50 - lam_x(k);
            50 + lam_x(k);
            50 - lam_y(k);
            50 + lam_y(k);
            ca(1, k)];
end
%%
% sol_eval = log.sol_approx;
sol_eval = log.sol;

ns = dof * (Ns+1);

offset = (dof + 1) * (Ns + 1);
x_num = sol_eval(1:offset);

y_num = sol_eval(offset+1:offset + ns);
offset = offset + ns;

c_num = sol_eval(offset+1:offset + ns);
offset = offset + ns;

s_num = sol_eval(offset+1:offset + ns);
offset = offset + ns;

vx_num = sol_eval(offset+1:offset + (dof + 1) * (Ns + 1));
offset = offset + (dof + 1) * (Ns + 1);

vy_num = sol_eval(offset+1:offset + ns);
offset = offset + ns;

ca_num = sol_eval(offset+1:offset + ns);
offset = offset + ns;

sa_num = sol_eval(offset+1:offset + ns);
offset = offset + ns;


lam_x_num = sol_eval(offset+1:offset+dof * Ns);
offset = offset + dof * Ns;

lam_y_num = sol_eval(offset+1:offset+dof * Ns);
offset = offset + dof * Ns;

tau_num = sol_eval(offset+1:offset+dof * Ns);
%%
tau_num(tau_num >  20) =  20;
tau_num(tau_num < -20) = -20;

lam_y_num(lam_y_num >  50) =  50;
lam_y_num(lam_y_num < -50) = -50;

lam_x_num(lam_x_num >  50) =  50;
lam_x_num(lam_x_num < -50) = -50;

norm_1 = sqrt(c_num.^2 + s_num.^2);
norm_2 = sqrt(ca_num.^2 + sa_num.^2);

c_num = c_num ./ norm_1;
s_num = s_num ./ norm_1;

ca_num = ca_num ./ norm_2;
sa_num = sa_num ./ norm_2;

assign(x, reshape(x_num, size(x)));
assign(y, reshape(y_num, size(y)));
assign(c, reshape(c_num, size(c)));
assign(s, reshape(s_num, size(s)));

assign(vx, reshape(vx_num, size(vx)));
assign(vy, reshape(vy_num, size(vy)));
assign(ca, reshape(ca_num, size(ca)));
assign(sa, reshape(sa_num, size(sa)));

assign(lam_x, reshape(lam_x_num, size(lam_x)));
assign(lam_y, reshape(lam_y_num, size(lam_y)));
assign(tau,   reshape(tau_num, size(tau)));
%%
% ops = sdpsettings('solver','fmincon');%,'fmincon.algorithm','sqp');
ops = sdpsettings('solver','ipopt');
ops.print_level = 12;
ops.verbose = true;

% ops.fmincon.MaxIterations = 1e7;
ops.ipopt.max_iter = 1e4 * 10;
ops.ipopt.max_cpu_time = 1e100;

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
sol_refine = [x_pole(:); x_cart(:); y_pole(:); c_pole(:); s_pole(:);...
              vx_pole(:); vx_cart(:); vy_pole(:); ca_pole(:); sa_pole(:);...
              lam_x(:); lam_y(:); tau(:)];
end