function [R_num, F_num, p_num, v_num, tau_num, f_num, p_num_load, v_num_load, lam_num, obj] = refine_solution_det_load(sol, dof, dt, mass, inertial, Qc, Rc, Pc, quat_init, obs_flag, dx0_in, dz0_in)
% constant
% Jb = diag([1, 3, 6]);
Jb = inertial;
m = mass;
m_load = 0.5;
arm_length = 0.5;

% dt = 0.0625;%25;
g = [0;0;-9.81];
R_init = quat2rotm(quat_init);

R_goal = eye(3);
%
yalmip('clear')
R = sdpvar(dof+1, 3, 3);
F = sdpvar(dof+1, 3, 3);

p = sdpvar(dof+1, 3, 1);
v = sdpvar(dof+1, 3, 1);

tau = sdpvar(dof, 3, 1);
f = sdpvar(dof, 1);

p_load = sdpvar(dof+1, 3, 1);
v_load = sdpvar(dof+1, 3, 1);

% load
lam    = sdpvar(dof, 1);


eq = [];
ineq = [];
f_act = 0;
% initial condition
dx0 = dx0_in; %2;
dz0 = dz0_in; % (arm_length - sqrt(( arm_length^2 - (dx0 * dt)^2 ))) / dt;
v0_load = [0; dx0; dz0];
eq = [eq;
    R(1,:,:) == R_init;
    F(1,:,:) == eye(3);
    p(1,:) == [1;1;3]';
    v(1,:) == [0;0; 0]';
    % load
    p_load(1,:) == [1;1;2.5]';
    v_load(1,:) == v0_load';
    ];

for k = 1:dof
    Rk  = squeeze(R(k, :, :));
    Rkp = squeeze(R(k+1,:,:));
    Fk  = squeeze(F(k, :, :));
    Fkp = squeeze(F(k+1,:,:));
    
    pk  = squeeze(p(k, :, :))';
    pkp = squeeze(p(k+1,:,:))';
    vk  = squeeze(v(k, :, :))';
    vkp = squeeze(v(k+1,:,:))';

    pk_load  = squeeze(p_load(k, :, :))';
    pkp_load = squeeze(p_load(k+1,:,:))';
    vk_load  = squeeze(v_load(k, :, :))';
    vkp_load = squeeze(v_load(k+1,:,:))';
    
    tauk = tau(k, :)';
    fk   = f(k);
    lamk = lam(k);
    % SO3 constraints
    IR = Rk' * Rk - eye(3);
    IF = Fk' * Fk - eye(3);
    eq = [eq;
        IR(1,1) == 0;
        IR(2,2) == 0;
        IR(3,3) == 0;
        IR(1,2) == 0;
        IR(1,3) == 0;
        IR(2,3) == 0;
        
        IF(1,1) == 0;
        IF(2,2) == 0;
        IF(3,3) == 0;
        IF(1,2) == 0;
        IF(1,3) == 0;
        IF(2,3) == 0];
    dR = Rkp - Rk * Fk;
    eq = [eq;
        dR(1,1) == 0;
        dR(2,2) == 0;
        dR(3,3) == 0;
        dR(1,2) == 0;
        dR(1,3) == 0;
        dR(2,3) == 0;
        %
        pkp - (pk + dt * Rk * vk) == 0;
        pkp_load - (pk_load + dt * vk_load) == 0];
    % dynamics constraints
    dM = Fkp * Jb - Jb * Fkp' - (Jb * Fk - Fk' * Jb);

    arm = pkp - pkp_load;

    dv      = m * vkp - (m * Fk' * vk + ( fk * [0;0;1] + Rkp' * (m * g + lamk * arm)  ) * dt );
    dv_load = m_load * vkp_load - (m_load * vk_load + (-lamk * arm + m_load * g) * dt );
    
    eq = [eq;
        dM(1, 2) + (-tauk(3)) * dt^2 == 0;
        dM(1, 3) + ( tauk(2)) * dt^2 == 0;
        dM(2, 3) + (-tauk(1)) * dt^2 == 0;
        dv == 0;
        dv_load == 0;
        sum(arm.^2) - arm_length.^2 == 0];
    
    p_goal = [0, 0, 0.5]';
    p_goal_load = [0, 0, 0]';

    R_cost = Qc(1) * trace((Rk - R_goal)' * (Rk - R_goal));
    w_cost = Qc(2) * trace( (Fk - eye(3))' * (Fk - eye(3)) );
    p_cost = Qc(3) * sum((pk - p_goal).^2);
    v_cost = Qc(4) * sum(vk.^2);
    
    p_load_cost = Qc(3) * sum((pk_load - p_goal_load).^2);
    v_load_cost = Qc(4) * sum(vk_load.^2);

    ctr_cost = Rc * sum(tauk.^2) + Rc* fk^2;
    f_act = f_act + R_cost + w_cost + p_cost + v_cost + ctr_cost + p_load_cost + v_load_cost;
    
    ineq = [ineq;
            -5 <= tauk;
            tauk <= 5;

            sum((pkp - pk).^2) <= 0.45^2;
            sum((pkp_load - pk_load).^2) <= 0.45^2;       
%             -7 <= fk;
%             fk <= 7;
%             (pkp(2) - 0.5)^2 + (pkp(1) - 0.5)^2 >= 0.25;
%             (pkp(1))^2 + (pkp(2) - 0.5)^2 >= 0.25;
%             (pkp(1) - 0.6)^2 + (pkp(2) - 0.5)^2 >= 0.25;
            pkp(3) >= 0];
    if obs_flag
        ineq = [ineq; 
                (pkp(1) - 0.6)^2 + (pkp(2) - 0.5)^2 >= 0.25;
                (pkp_load(1) - 0.6)^2 + (pkp_load(2) - 0.5)^2 >= 0.25];
    end
end

R_cost_term = Pc(1) * trace((Rkp - R_goal)' * (Rkp - R_goal));
F_cost_term = Pc(2) * trace((Fkp - eye(3))' * (Fkp - eye(3)));
p_cost_term = Pc(3) * sum((pkp - p_goal).^2);
v_cost_term = Pc(4) * sum(vkp.^2);

p_cost_term_load = Pc(3) * sum((pkp_load - p_goal_load).^2);
v_cost_term_load = Pc(4) * sum(vkp_load.^2);


cost = f_act * dt / 0.25 + R_cost_term + F_cost_term + p_cost_term + v_cost_term;
cost = cost + p_cost_term_load + v_cost_term_load;
cost = cost * 10;
%%
dof = double(dof);
R_num = sol(1:9 * (dof + 1));
F_num = sol(1+9 * (dof + 1):9 * (dof + 1) * 2);

shift = 9 * (dof + 1) * 2;
p_num = sol(          1+shift : shift+3*(dof+1));
v_num = sol(1+shift+3*(dof+1) : shift+6*(dof+1));
p_num = reshape(p_num, [dof+1, 3]);
v_num = reshape(v_num, [dof+1, 3]);

shift = shift+6*(dof+1);
control = sol(shift+1:shift + 4 * dof);

shift   = shift+4*dof;
p_num_load   = sol(          1+shift : shift+3*(dof+1));
v_num_load   = sol(1+shift+3*(dof+1) : shift+6*(dof+1));
p_num_load = reshape(p_num_load, [dof+1, 3]);
v_num_load = reshape(v_num_load, [dof+1, 3]);

shift = shift+6*(dof+1);
lam_num = sol(shift+1:end);

R_num = reshape(R_num, [dof+1, 3, 3]);
F_num = reshape(F_num, [dof+1, 3, 3]);  

%%% project to SO(3) %%%
% for k = 1:dof+1
%     Rk_before = squeeze(R_num(k,:,:));
%     [U, ~, V] = svd(Rk_before);
%     Rk_hat = U * diag([1, 1, det(U) * det(V)]) * V';
%     R_num(k,:,:) = Rk_hat;
% 
%     Fk_before = squeeze(F_num(k,:,:));
%     [U, ~, V] = svd(Fk_before);
%     Fk_hat = U * diag([1, 1, det(U) * det(V)]) * V';
%     F_num(k,:,:) = Fk_hat;
% end
%%% project to SO(3) %%% 


tau_num = reshape(control(1:dof*3), [dof, 3]);
tau_num(tau_num >  5) =  5;
tau_num(tau_num < -5) = -5;
f_num = reshape(control(dof*3+1:end), [dof, 1]);
f_num(f_num >  30) =  30;
f_num(f_num < -30) = -30;

assign(R, R_num)
assign(F, F_num)
assign(p, p_num)
assign(v, v_num)
assign(tau, tau_num)
assign(f, f_num)

assign(p_load, p_num_load)
assign(v_load, v_num_load)
assign(lam, lam_num)

%%
% ops = sdpsettings('solver','fmincon');%,'fmincon.algorithm','sqp');
ops = sdpsettings('solver','ipopt');
ops.print_level = 12;
ops.verbose = true;

% ops.fmincon.MaxIterations = 1e7;
ops.ipopt.max_iter = 1e4 * 10;
ops.ipopt.max_cpu_time = 1e100;

ops.usex0 = 1;
optimize([eq; ineq],cost,ops);
%%
R_num = value(R);
F_num = value(F);
p_num = value(p);
v_num = value(v);
tau_num = value(tau);

p_num_load = value(p_load);
v_num_load = value(v_load);
lam_num = value(lam);

f_num = value(f);
obj = value(cost);
end