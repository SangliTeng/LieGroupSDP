module Control

mutable struct param_base
    dt::Float64
    Ns::Int
    order::Int
#     CS_type::String
#     TS_type::String
end

mutable struct param_cart_multi
    control_param::param_base
    mass::Float64
    inertial::Float64
    lx::Float64
    ly::Float64
    run_cost::Matrix{Float64}
    terminal_cost::Matrix{Float64}
    terminal_cons::Bool #
    # x2 x1 y1 (c1 s1) || vx2 vx1 vy1 (ca1 sa1) || tau
    xc_0::Float64
    theta_0::Matrix{Float64}
    dof::Int
end


mutable struct param_cartpole
    control_param::param_base
    mass::Float64
    inertial::Float64
    lx::Float64
    ly::Float64
    run_cost::Matrix{Float64}
    terminal_cost::Matrix{Float64}
    terminal_cons::Bool #
    # x2 x1 y1 (c1 s1) || vx2 vx1 vy1 (ca1 sa1) || tau
    xc_0::Float64
    theta_0::Float64
end


mutable struct param_drone_landing
    control_param::param_base
    obs_num::Int
    obs_pos::Matrix{Float64}
end

end