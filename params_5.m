%load parameters
run aerosonde_parameters.m

%computition of trim conditions 
Va=30; %flight speed
R=9999999; %turning radius
gamma=15*pi/180; %flight path angle


%computition of trim point
[x_trim,u_trim,y_trim,dx_trim] = compute_trim('mavsim_trim',Va,R,gamma);

MAV.pe0=    0;
MAV.pn0=    0;
MAV.pd0=    0;
MAV.u0=     x_trim(4);
MAV.v0=     x_trim(5);
MAV.w0=     x_trim(6);
MAV.phi0   = x_trim(7);     % initial roll angle
MAV.theta0 = -90*pi/180;     % initial pitch angle
MAV.psi0   = x_trim(9);     % initial yaw angle
MAV.p0     = x_trim(10);     % initial body frame roll rate
MAV.q0     = x_trim(11);     % initial body frame pitch rate
MAV.r0     = x_trim(12);
MAV.delta_e0= u_trim(1);
MAV.delta_a0= u_trim(2);
MAV.delta_r0= u_trim(3);
MAV.delta_t0= u_trim(4);

%computition of transfer functions

[T_phi_delta_a, T_chi_phi, T_theta_delta_e, T_h_theta, T_h_Va,...
    T_Va_delta_t, T_Va_theta, T_v_delta_r] = compute_tf_model('mavsim_trim',x_trim,u_trim,MAV);


%computition of state space system

[A_lon, B_lon, A_lat, B_lat]=compute_ss_model('mavsim_trim',x_trim,u_trim);


