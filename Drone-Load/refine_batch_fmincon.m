clear all;clc;close all;
%%
for k = 5

    "k = " + num2str(k)
    file_name = "log_NC_" + num2str(k) + ".mat";
    log = load(file_name);

    dof = double(log.dof);
    dx0 = double(log.dx0);
    dz0 = double(log.dz0);

    sol = log.sol;
    dt = log.dt;
    Mass = log.Mass;
    Inertial = log.Inertial;
    Qc = log.Qc; 
    Rc = log.Rc;
    Pc = log.Pc;
    quat_init = log.quat_init; 
    obs_flag = log.obs_flag;

    [R, F, p, v, tau, f, p_load, v_load, lam, opt_fmincon] = refine_solution_det_load(sol + randn(size(sol)) * 0e-2, dof, dt, Mass, Inertial, Qc, Rc, Pc, double(quat_init'), obs_flag, dx0, dz0);
    "k = " + num2str(k) + " completed"
    file_name = "refine_NC_" + num2str(k) + ".mat";
    save_file(file_name, R, F, p, v, tau, f, p_load, v_load, lam, opt_fmincon)
end
%%
function save_file(file_name, R, F, p, v, tau, f, p_load, v_load, lam, opt_fmincon)
    save(file_name, "R", "F", "p", "v", "tau", "f", "p_load", "v_load", "lam", "opt_fmincon")
end