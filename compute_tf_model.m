% x_trim is the trimmed state,
% u_trim is the trimmed input
function [T_phi_delta_a, T_chi_phi, T_theta_delta_e, T_h_theta, T_h_Va,...
    T_Va_delta_t, T_Va_theta, T_v_delta_r] = compute_tf_model(filename,x_trim,u_trim,MAV);
Va_trim=x_trim(4);
Va=MAV.Va0;
theta_trim=x_trim(8);
%% paramaters
a_phi1=-0.5*MAV.rho*Va^2*MAV.S_wing*MAV.b*MAV.C_p_p*MAV.b/(2*Va);
a_phi2=0.5*MAV.rho*Va^2*MAV.S_wing*MAV.b*MAV.C_p_delta_a;

a_beta1=-MAV.rho*Va*MAV.S_wing/(2*MAV.mass)*MAV.C_Y_beta;
a_beta2=MAV.rho*Va*MAV.S_wing/(2*MAV.mass)*MAV.C_Y_delta_r;

a_theta1=-MAV.rho*Va*MAV.S_wing*MAV.c/(2*MAV.Jy)*MAV.C_m_q*MAV.c/(2*MAV.mass);
a_theta2=-MAV.rho*Va*MAV.S_wing*MAV.c/(2*MAV.Jy)*MAV.C_m_alpha;
a_theta3=MAV.rho*Va^2*MAV.c*MAV.S_wing/(2*MAV.Jy)*MAV.C_m_delta_e;

a_V1=MAV.rho*MAV.u0*MAV.S_wing/MAV.mass*(MAV.C_D_0+MAV.C_D_alpha*MAV.alpha0+...
    MAV.C_D_delta_e*MAV.delta_e0)+MAV.rho*MAV.S_prop*MAV.C_prop*MAV.u0/MAV.mass;
a_V2=MAV.rho*MAV.S_prop*MAV.C_prop*MAV.k_motor*MAV.delta_t0;
a_V3=MAV.gravity*cos(MAV.theta0-MAV.psi0);
%%
% define transfer functions
T_phi_delta_a   = tf([a_phi2],[1,a_phi1,0]);
T_chi_phi       = tf([MAV.gravity/Va_trim],[1,0]);
T_theta_delta_e = tf(a_theta3,[1,a_theta1,a_theta2]);
T_h_theta       = tf([Va_trim],[1,0]);
T_h_Va          = tf([theta_trim],[1,0]);
T_Va_delta_t    = tf([a_V2],[1,a_V1]);
T_Va_theta      = tf([-a_V3],[1,a_V1]);
T_v_delta_r     = tf([a_beta2],[1,a_beta1]);
