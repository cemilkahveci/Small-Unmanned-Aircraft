% forces_moments.m
%   Computes the forces and moments acting on the airframe. 
%
%   Output is
%       F     - forces
%       M     - moments
%       Va    - airspeed
%       alpha - angle of attack
%       beta  - sideslip angle
%       wind  - wind vector in the inertial frame
%

function out = forces_moments(x, delta, wind,MAV)

    % relabel the inputs
    pn      = x(1);
    pe      = x(2);
    pd      = x(3);
    u       = x(4);
    v       = x(5);
    w       = x(6);
    phi     = x(7);
    theta   = x(8);
    psi     = x(9);
    p       = x(10);
    q       = x(11);
    r       = x(12);
    delta_e = delta(1);
    delta_a = delta(2);
    delta_r = delta(3);
    delta_t = delta(4);
    w_ns    = wind(1); % steady wind - North
    w_es    = wind(2); % steady wind - East
    w_ds    = wind(3); % steady wind - Down
    u_wg    = wind(4); % gust along body x-axis
    v_wg    = wind(5); % gust along body y-axis    
    w_wg    = wind(6); % gust along body z-axis
    
    % compute wind data in NED
    w_gs=invrotate([u_wg;v_wg;w_wg],phi,theta,psi);
    
    w_n = w_gs(1)+w_ns;
    w_e = w_gs(2)+w_es;
    w_d = w_gs(3)+w_ds;
    
    vbw=rotate([w_ns;w_es;w_ds],phi,theta,psi)+[u_wg;v_wg;w_wg];
    
    ur=u-vbw(1);
    vr=v-vbw(2);
    wr=w-vbw(3);

    % compute air data
    Va = sqrt(ur^2+vr^2+wr^2);
    alpha = atan(wr/ur);
    beta = asin(vr/sqrt(ur^2+vr^2+wr^2));
    
%     % compute external forces and torques on aircraft
%     transformation_matrix=[cos(alpha) -sin(alpha);sin(alpha) cos(alpha)];
     f_drag=0.5*MAV.rho*Va^2*MAV.S_wing*(MAV.C_D_0+MAV.C_D_alpha*alpha+...
         MAV.C_D_q*MAV.c/(2*Va)*q+MAV.C_D_delta_e*delta_e);
     f_lift=0.5*MAV.rho*Va^2*MAV.S_wing*(MAV.C_L_0+MAV.C_L_alpha*alpha+...
         MAV.C_L_q*MAV.c/(2*Va)*q+MAV.C_L_delta_e*delta_e);

    Force(1) =-f_drag*cos(alpha)+f_lift*sin(alpha) +...
        0.5*MAV.rho*MAV.S_prop*MAV.C_prop*((MAV.k_motor*delta_t)^2-Va^2)-...
        MAV.mass*MAV.gravity*sin(theta);
    Force(2) =  0.5*MAV.rho*Va^2*MAV.S_wing*(MAV.C_Y_0+MAV.C_Y_beta*beta+...
        MAV.C_Y_p*MAV.c/(2*Va)*p+MAV.C_Y_r*MAV.c/(2*Va)*r+MAV.C_Y_delta_a*delta_a+...
        MAV.C_Y_delta_r*delta_r)+...
        MAV.mass*MAV.gravity*cos(theta)*sin(phi);
    Force(3) =-f_drag*sin(alpha)-f_lift*cos(alpha) +...
        MAV.mass*MAV.gravity*cos(theta)*cos(phi);

    Torque(1) = 0.5*MAV.rho*Va^2*MAV.S_wing*MAV.b*(MAV.C_ell_0+MAV.C_ell_beta*beta+...
        MAV.C_ell_p*MAV.b/(2*Va)*p+MAV.C_ell_r*MAV.b/(2*Va)*r+MAV.C_ell_delta_a*delta_a+...
        MAV.C_ell_delta_r*delta_r)-MAV.k_T_P*(MAV.k_Omega*delta_t)^2;
    Torque(2) = 0.5*MAV.rho*Va^2*MAV.S_wing*MAV.c*(MAV.C_m_0+MAV.C_m_alpha*alpha+...
        MAV.C_m_q*MAV.c/(2*Va)*q+MAV.C_m_delta_e*delta_e);
    Torque(3) = 0.5*MAV.rho*Va^2*MAV.S_wing*MAV.b*(MAV.C_n_0+MAV.C_n_beta*beta+...
        MAV.C_n_p*MAV.b/(2*Va)*p+MAV.C_n_r*MAV.b/(2*Va)*r+MAV.C_n_delta_a*delta_a+...
        MAV.C_n_delta_r*delta_r);
   
    out = [Force'; Torque'; Va; alpha; beta; w_n; w_e; w_d];
end


