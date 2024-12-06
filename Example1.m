%% فرضیات و داده های سوال
diameter = 0.0762; % m
q_oil = 9.380; % m^3/d
q_oil = q_oil/86400; % m3/s
q_water = 86.17; 
q_water =q_water/86400;
q_gas = 1161; 
q_gas = q_gas/86400;
den_ref_oil = 853.4; % kg/m3
den_ref_water = 1010;  
den_ref_gas = 1.32; 
gravity = 10;
teta = 0; % 
dt = 100; % 100 s
dz = 200; % m
C0 = 1;
V_d = 0;
Area = 2*pi*(diameter/2)*dz;
volume = pi*(diameter/2)^2*dz;
Area_V = pi*(diameter/2)^2;
T_top = 24.44 + 273.15;
T_bottom = 42.22 + 273.15;
depth = 1632;
number_grids = ceil (depth/dz) ;
U_to_top = 100;  % W/m2 C
U_to_end = 110;
Cp_W = 4.18; 
Cp_O = 2; 
Cp_G = 1; 

%% density and residual handle functions
rho_P = @ (p_ref, p, compressibility, den_ref) den_ref *exp(compressibility* (p-p_ref));
rho_T = @ (T_ref, T, thermal_expansion, den_ref) den_ref *(1+thermal_expansion* (T-T_ref));
rho_TP = @ (rho_P, rho_T, den_ref) rho_P* rho_T/den_ref;
mixure_density = @ (alpha_g, alpha_w, density_gas, density_water, density_oil) ...
    (alpha_g*density_gas+ alpha_w*density_water + (1-alpha_g - alpha_w)*density_oil);

% ***************** Momentum Equation *****************
Residiual1 = @ (pressure_i_before, pressure_i_present, mixure_density_next_x_time,...
    mixure_density_present, mixure_density_next_time, mixture_velocity_next_x_time, ...
    mixture_velocity_present, mixture_velocity_next_time, f_tp) (pressure_i_present - ...
    pressure_i_before) + mixure_density_next_time*gravity*cosd (teta) + (mixure_density_next_time *...
    mixture_velocity_next_time -  mixure_density_present * mixture_velocity_present) / dt + ...
    (mixure_density_next_x_time * mixture_velocity_next_x_time^2 - mixure_density_next_time *...
    mixture_velocity_next_time^2) / dz + (f_tp * mixure_density_next_time * mixture_velocity_next_time^2)/(2*diameter);

% ***************** Mass conservation (Water) *****************
Residiual2 = @ (water_density_present, water_density_next_time, water_density_next_x_time, ...
    alpha_water_present, alpha_water_next_time, alpha_water_next_x_time, V_w, q_ww) ...
    (dz*Area/dt)*(water_density_next_time* alpha_water_next_time - water_density_present* ...
    alpha_water_present) - Area*(water_density_next_x_time*alpha_water_next_x_time*V_w - ...
    water_density_next_time*alpha_water_next_time*V_w) - water_density_next_time*q_ww;

% ***************** Mass conservation (Oil) *****************
Residiual3 = @ (gas_density_present, gas_density_next_time, gas_density_next_x_time, ...
    alpha_gas_present, alpha_gas_next_time, alpha_gas_next_x_time, oil_density_present, ...
    oil_density_next_time, oil_density_next_x_time, alpha_oil_present, alpha_oil_next_time,...
    alpha_oil_next_x_time, Xog, Xoo, qoo, qog, gas_velocity_next_x_time, gas_velocity_next_time,...
    oil_velocity_next_x_time, oil_velocity_next_time) (dz*Area/dt)*(gas_density_next_time*...
    alpha_gas_next_time*Xog + oil_density_next_time*alpha_oil_next_time*Xoo - ...
    gas_density_present*alpha_gas_present*Xog - oil_density_present*alpha_oil_present*Xoo) ...
    - Area*(gas_density_next_x_time*alpha_gas_next_x_time*gas_velocity_next_x_time*Xog +...
    oil_density_next_x_time* alpha_oil_next_x_time*oil_velocity_next_x_time*Xoo - ...
    - gas_density_next_time*alpha_gas_next_time*gas_velocity_next_time*Xog -...
    oil_density_next_time*alpha_oil_next_time*oil_velocity_next_time*Xoo) - ...
    oil_density_next_time*qoo - gas_density_next_time*qog;

% ***************** Mass conservation (Gas) *****************
Residiual4 = @ (gas_density_present, gas_density_next_time, gas_density_next_x_time, ...
    alpha_gas_present, alpha_gas_next_time, alpha_gas_next_x_time, oil_density_present,...
    oil_density_next_time, oil_density_next_x_time, alpha_oil_present,alpha_oil_next_time,...
    alpha_oil_next_x_time,water_density_present, water_density_next_time, ...
    water_density_next_x_time, alpha_water_present, alpha_water_next_time, ...
    alpha_water_next_x_time, Xgg, Xgo, Xgw, qgg, qgo, qgw,gas_velocity_next_x_time,...
    gas_velocity_next_time, oil_velocity_next_x_time, oil_velocity_next_time,...
    water_velocity_next_x_time, water_velocity_next_time) (dz*Area/dt)*(gas_density_next_time*...
    alpha_gas_next_time*Xgg + oil_density_next_time*alpha_oil_next_time*Xgo + ...
    water_density_next_time*alpha_water_next_time*Xgw - gas_density_present*alpha_gas_present*Xgg -...
    oil_density_present*alpha_oil_present*Xgo - water_density_present*alpha_water_present*Xgw) ...
    - Area*(gas_density_next_x_time*alpha_gas_next_x_time*gas_velocity_next_x_time*Xgg +...
    oil_density_next_x_time*alpha_oil_next_x_time*oil_velocity_next_x_time*Xgo - ...
    water_density_next_x_time*alpha_water_next_x_time*water_velocity_next_x_time*Xgw ...
    - gas_density_next_time*alpha_gas_next_time*gas_velocity_next_time*Xgg -...
    oil_density_next_time*alpha_oil_next_time*oil_velocity_next_time*Xgo - ...
    water_density_next_time*alpha_water_next_time*water_velocity_next_time*Xgo) - ...
    oil_density_next_time*qgo - gas_density_next_time*qgg  - water_density_next_time*qgw;

% ***************** Energy Balance *****************
Residiual5 = @ (gas_density_present, gas_density_next_time, gas_density_next_x_time, ...
    alpha_gas_present, alpha_gas_next_time, alpha_gas_next_x_time, oil_density_present, ...
    oil_density_next_time, oil_density_next_x_time, alpha_oil_present, alpha_oil_next_time, ...
    alpha_oil_next_x_time, water_density_present, water_density_next_time, ...
    water_density_next_x_time, alpha_water_present, alpha_water_next_time, ...
    alpha_water_next_x_time, gas_velocity_next_x_time, gas_velocity_next_time, ...
    gas_velocity_present, oil_velocity_next_x_time, oil_velocity_next_time, oil_velocity_present, ...
    water_velocity_next_x_time, water_velocity_next_time, water_velocity_present , ...
    Upg_present, Upg_next_time, Upg_next_x_time, Upo_present, Upo_next_time, ...
    Upo_next_x_time, Upw_present, Upw_next_time, Upw_next_x_time ,pressure_next_time, ...
    pressure_next_x_time,Q_loss) (dz*Area/dt)*(gas_density_next_time*alpha_gas_next_time* ...
    (Upg_next_time+0.5* gas_velocity_next_time)+ oil_density_next_time*alpha_oil_next_time*...
    (Upo_next_time+0.5* oil_velocity_next_time) + water_density_next_time*alpha_water_next_time*...
    (Upw_next_time+0.5* water_velocity_next_time) - gas_density_present*alpha_gas_present*...
    (Upg_present+0.5* gas_velocity_present) - oil_density_present*alpha_oil_present* ...
    (Upo_present+0.5* oil_velocity_present) - water_density_present*alpha_water_present*...
    (Upw_present+0.5* water_velocity_present)) - Area*(gas_density_next_x_time* ...
    alpha_gas_next_x_time*gas_velocity_next_x_time*(Upg_next_x_time+pressure_next_x_time/ gas_density_next_x_time)...
    +oil_density_next_x_time*alpha_oil_next_x_time*oil_velocity_next_x_time*(Upo_next_x_time+...
    pressure_next_x_time/ oil_density_next_x_time) + water_density_next_x_time* ...
    alpha_water_next_x_time*water_velocity_next_x_time* (Upw_next_x_time+ ...
    pressure_next_x_time/ water_density_next_x_time) - gas_density_next_time* ...
    alpha_gas_next_time*gas_velocity_next_time*(Upg_next_time+pressure_next_time/ gas_density_next_time)...
    -oil_density_next_time*alpha_oil_next_time*oil_velocity_next_time*(Upo_next_time+pressure_next_time/ oil_density_next_time)...
     - water_density_next_time*alpha_water_next_time*water_velocity_next_time*(Upw_next_time+ ...
     pressure_next_time/ water_density_next_time)) - volume*(gas_density_next_time* ...
     alpha_gas_next_time*gas_velocity_next_time*gravity*cosd(teta) + oil_density_next_time*...
     alpha_oil_next_time*oil_velocity_next_time*gravity*cosd(teta) + water_density_next_time*...
     alpha_water_next_time*water_velocity_next_time*gravity*cosd(teta)) + dz*Q_loss;
    

iteration_number = 1;
for time = 0: dt: 864*dt % 1 days
    if time == 0 && iteration_number ==1   % initial condition
        for i = 1: number_grids   
            initial_T (i,1) = reservoir_temp (T_top, T_bottom, i) + 4; %#ok
            initial_P (i,1) = (0.6*101325/(0.3048*14.7)) * i * dz; %#ok
            initial_alpha_gas (i,1) = 0.3; %#ok
            initial_alpha_water (i,1) = 0.2; %#ok
            initial_alpha_oil (i,1) = 1-initial_alpha_gas (i,1)-initial_alpha_water(i,1); %#ok
            initial_V_m (i,1) = initial_alpha_gas (i,1).*(q_gas/Area_V) + initial_alpha_water (i,1).*(q_water/Area_V) + initial_alpha_oil (i,1).*(q_oil/Area_V); %#ok
        end
        delta_T1 = 2;
        delta_P1 = +10;
        delta_alpha_g1 = -0.1;
        delta_alpha_w1 = -0.1;
        delta_alpha_Vm1 = 0.01;
        delta1 = [delta_T1;delta_P1;delta_alpha_g1;  delta_alpha_w1; delta_alpha_Vm1];
        initial_guess_T =  initial_T (:,1) + delta_T1;
        initial_guess_P = initial_P (:,1) + delta_P1;
        initial_guess_alpha_g = initial_alpha_gas (:,1) + delta_alpha_g1;
        initial_guess_alpha_w = initial_alpha_water (:,1) + delta_alpha_w1;
        initial_guess_alpha_o = 1-  initial_guess_alpha_g(:,1)- initial_guess_alpha_w(:,1);
        initial_guess_V_m = initial_V_m (:,1) + delta_alpha_Vm1;
         initial_V_m (1,1) = 0;
    end


    while true        
        for j = 1: number_grids
            if j == 1 
                initial_P (1,1)=779100;
                initial_guess_P (1,1) = 779100;
            elseif j == number_grids
                initial_guess_T (j+1,1) = reservoir_temp (T_top, T_bottom, j);
                initial_guess_P (j+1,1) = initial_guess_P (j,1);
                initial_guess_alpha_g (j+1,1) = initial_guess_alpha_g (j,1);
                initial_guess_alpha_w (j+1,1) =  initial_guess_alpha_w (j,1);
                initial_guess_alpha_o (j+1,1) =  initial_guess_alpha_o (j,1);
                initial_guess_V_m(j+1,1) = initial_guess_V_m(j,1);
            end
                
            % density of water i n
            % p ref and p in psi. c-1 is psi-1.
            rho_water_in = rho_P (14.7, (14.7*initial_P (j,1)/101325), 10^-6, den_ref_water) * rho_T (20+273.15,  initial_T (j,1), 210*10^-6, den_ref_water) / den_ref_water;
            rho_water_in1 = rho_P (14.7, (14.7*initial_guess_P (j,1)/101325), 10^-6, den_ref_water) * rho_T (20+273.15, initial_guess_T (j,1), 210*10^-6, den_ref_water) / den_ref_water;
            rho_water_i1n1 = rho_P (14.7, (14.7*initial_guess_P (j+1,1)/101325), 10^-6, den_ref_water) * rho_T (20+273.15, initial_guess_T (j+1,1), 210*10^-6, den_ref_water) / den_ref_water;
            rho_oil_in = rho_P (14.7, (14.7*initial_P (j,1)/101325), 6*10^-6, den_ref_oil) * rho_T (20+273.15,  initial_T (j,1), 950*10^-6, den_ref_oil) / den_ref_oil;
            rho_oil_in1 = rho_P (14.7, (14.7*initial_guess_P (j,1)/101325), 6*10^-6, den_ref_oil) * rho_T (20+273.15, initial_guess_T (j,1), 950*10^-6, den_ref_oil) / den_ref_oil;
            rho_oil_i1n1 = rho_P (14.7, (14.7*initial_guess_P (j+1,1)/101325), 6*10^-6, den_ref_oil) * rho_T (20+273.15, initial_guess_T (j+1,1), 950*10^-6, den_ref_oil) / den_ref_oil;                  
            rho_gas_in = rho_P (14.7, (14.7*initial_P (j,1)/101325), 100*10^-6, den_ref_gas) * rho_T (20+273.15,  initial_T (j,1), 3400*10^-6, den_ref_gas) / den_ref_gas;
            rho_gas_in1 = rho_P (14.7, (14.7*initial_guess_P (j,1)/101325), 100*10^-6, den_ref_gas) * rho_T (20+273.15, initial_guess_T (j,1), 3400*10^-6, den_ref_gas) / den_ref_gas;
            rho_gas_i1n1 = rho_P (14.7, (14.7*initial_guess_P (j+1,1)/101325), 100*10^-6, den_ref_gas) * rho_T (20+273.15, initial_guess_T (j+1,1), 3400*10^-6, den_ref_gas) / den_ref_gas;         
            rho_mixture_in = mixure_density (initial_alpha_gas (j,1), initial_alpha_water (j,1), rho_gas_in, rho_water_in,  rho_oil_in);
            rho_mixture_in1 = mixure_density (initial_guess_alpha_g (j,1), initial_guess_alpha_w (j,1), rho_gas_in1, rho_water_in1,  rho_oil_in1);
            rho_mixture_i1n1 = mixure_density (initial_guess_alpha_g (j+1,1), initial_guess_alpha_w (j+1,1),  rho_gas_i1n1, rho_water_i1n1,  rho_oil_i1n1);
            f_tp = friction_tp (initial_guess_alpha_g (j,1), initial_guess_alpha_w (j,1), rho_mixture_in1, initial_guess_V_m (j,1));
            if j ~= 1
                R1_base = Residiual1 (initial_guess_P (j-1,1), initial_guess_P (j,1), rho_mixture_i1n1, rho_mixture_in, rho_mixture_in1, initial_guess_V_m (j+1,1) , initial_V_m (j,1), initial_guess_V_m (j,1), f_tp) ;
            else
                R1_base = Residiual1 (779100, initial_guess_P (j,1), rho_mixture_i1n1, rho_mixture_in, rho_mixture_in1, initial_guess_V_m (j+1,1) , initial_V_m (j,1), initial_guess_V_m (j,1), f_tp) ;
            end
            
            R2_base = Residiual2 (rho_water_in, rho_water_in1, rho_water_i1n1, initial_alpha_water (j,1), initial_guess_alpha_w (j,1), initial_guess_alpha_w (j+1,1), q_water/Area_V, q_water);            
            gas_velocity_present = gas_velocity (C0, V_d, initial_V_m(j,1));
            gas_velocity_next_time = gas_velocity (C0, V_d, initial_guess_V_m(j,1));
            gas_velocity_next_x_time = gas_velocity (C0, V_d, initial_guess_V_m(j+1,1));

            R3_base = Residiual3 (rho_gas_in, rho_gas_in1, rho_gas_i1n1, initial_alpha_gas (j,1), initial_guess_alpha_g(j,1), initial_guess_alpha_g (j+1,1), ...
            rho_oil_in, rho_oil_in1, rho_oil_i1n1,  initial_alpha_oil (j,1), initial_guess_alpha_o (j,1), initial_guess_alpha_o (j+1,1), 0, 1, q_oil, 0,...
            gas_velocity_next_x_time, gas_velocity_next_time, initial_guess_V_m(j+1,1)-gas_velocity_next_x_time-q_water/Area_V, initial_guess_V_m(j,1)-gas_velocity_next_time-q_water/Area_V);
        
            R4_base = Residiual4 (rho_gas_in, rho_gas_in1, rho_gas_i1n1, initial_alpha_gas (j,1), initial_guess_alpha_g (j,1), initial_guess_alpha_g (j+1,1), ...
            rho_oil_in, rho_oil_in1, rho_oil_i1n1, initial_alpha_oil (j,1), initial_guess_alpha_o (j,1), initial_guess_alpha_o (j+1,1),...
            rho_water_in, rho_water_in1, rho_water_i1n1, initial_alpha_water (j,1), initial_guess_alpha_w (j,1), initial_guess_alpha_w (j+1,1), 1, 0, 0, q_gas, 0, 0,...
            gas_velocity_next_x_time, gas_velocity_next_time, initial_guess_V_m(j+1,1)-gas_velocity_next_x_time-q_water/Area_V, initial_guess_V_m(j,1)-gas_velocity_next_time-q_water/Area_V, q_water/Area_V, q_water/Area_V) ;
             
            R5_base = Residiual5 (rho_gas_in, rho_gas_in1, rho_gas_i1n1, initial_alpha_gas (j,1), initial_guess_alpha_g (j,1), initial_guess_alpha_g (j+1,1), ...
                rho_oil_in, rho_oil_in1, rho_oil_i1n1, initial_alpha_oil (j,1), initial_guess_alpha_o (j,1), initial_guess_alpha_o (j+1,1),...
                rho_water_in, rho_water_in1, rho_water_i1n1, initial_alpha_water (j,1), initial_guess_alpha_w (j,1), initial_guess_alpha_w (j+1,1),...
                gas_velocity_next_x_time, gas_velocity_next_time, gas_velocity_present, initial_guess_V_m(j+1,1)-gas_velocity_next_x_time-q_water/Area, initial_guess_V_m(j,1)-gas_velocity_next_time-q_water/Area_V, initial_guess_V_m(j,1)-gas_velocity_present-q_water/Area_V,q_water/Area_V, q_water/Area_V, q_water/Area_V, ...
                Up (Cp_G, initial_T (j,1), 20+273.15), Up (Cp_G, initial_guess_T (j,1), 20+273.15), Up (Cp_G, initial_guess_T (j+1,1), 20+273.15),...
                Up (Cp_O, initial_T (j,1), 20+273.15), Up (Cp_O, initial_guess_T (j,1), 20+273.15), Up (Cp_O, initial_guess_T (j+1,1), 20+273.15),...
                Up (Cp_W, initial_T (j,1), 20+273.15), Up (Cp_W, initial_guess_T (j,1), 20+273.15), Up (Cp_W, initial_guess_T (j+1,1), 20+273.15),...
                initial_guess_P (j,1), initial_guess_P (j+1,1), Q_loss (T_top, T_bottom, U_to_top, U_to_end, initial_guess_T (j,1), j));
                Residual_mat = [R1_base; R2_base;R3_base; R4_base; R5_base];
                cell2{j}=Residual_mat;

        
            for k1 = 1: 5 % residual
                for k2 = 1:5 % derivative with respect to
                    if k1 == 1 && k2 == 1        
                        epsilon = 10^-2;
                        rho_mixture_in_1 = mixure_density (initial_alpha_gas (j,1)+epsilon, initial_alpha_water (j,1), rho_gas_in, rho_water_in,  rho_oil_in);
                        rho_mixture_in1_1 = mixure_density (initial_guess_alpha_g (j,1)+epsilon, initial_guess_alpha_w (j,1), rho_gas_in1, rho_water_in1,  rho_oil_in1);
                        rho_mixture_i1n1_1 = mixure_density (initial_guess_alpha_g (j+1,1), initial_guess_alpha_w (j+1,1),  rho_gas_i1n1, rho_water_i1n1,  rho_oil_i1n1);
                        f_tp_1 = friction_tp (initial_guess_alpha_g (j,1)+epsilon, initial_guess_alpha_w (j,1), rho_mixture_in1_1, initial_guess_V_m (j,1));
                        if j~=1
                            R1_1 = Residiual1 (initial_guess_P (j-1,1), initial_guess_P (j,1), rho_mixture_i1n1_1, rho_mixture_in_1, rho_mixture_in1_1, initial_guess_V_m (j+1,1) , initial_V_m (j,1), initial_guess_V_m (j,1), f_tp_1 ) ;
                        else
                            R1_1 = Residiual1 (779100, initial_guess_P (j,1), rho_mixture_i1n1_1, rho_mixture_in_1, rho_mixture_in1_1, initial_guess_V_m (j+1,1) , initial_V_m (j,1), initial_guess_V_m (j,1), f_tp_1 ) ;
                        end    
                        mat (k1, k2) = (R1_1 - R1_base)/epsilon;

                    elseif k1 == 1 && k2 == 2 
                        epsilon = 10^-2;
                        rho_mixture_in_2 = mixure_density (initial_alpha_gas (j,1), initial_alpha_water (j,1) + epsilon, rho_gas_in, rho_water_in,  rho_oil_in);
                        rho_mixture_in1_2 = mixure_density (initial_guess_alpha_g (j,1), initial_guess_alpha_w (j,1) + epsilon, rho_gas_in1, rho_water_in1,  rho_oil_in1);
                        rho_mixture_i1n1_2 = mixure_density (initial_guess_alpha_g (j+1,1), initial_guess_alpha_w (j+1,1),  rho_gas_i1n1, rho_water_i1n1,  rho_oil_i1n1);
                        f_tp_2 = friction_tp (initial_guess_alpha_g (j,1), initial_guess_alpha_w (j,1)+epsilon, rho_mixture_in1_2, initial_guess_V_m (j,1));
                        if j~=1
                            R1_2 = Residiual1 (initial_guess_P (j-1,1), initial_guess_P (j,1), rho_mixture_i1n1_2, rho_mixture_in_2, rho_mixture_in1_2, initial_guess_V_m (j+1,1) , initial_V_m (j,1), initial_guess_V_m (j,1), f_tp_2 ) ;
                        else
                            R1_2 = Residiual1 (779100, initial_guess_P (j,1), rho_mixture_i1n1_2, rho_mixture_in_2, rho_mixture_in1_2, initial_guess_V_m (j+1,1) , initial_V_m (j,1), initial_guess_V_m (j,1), f_tp_2 ) ;
                        end               
                        mat (k1, k2) = (R1_2 - R1_base)/epsilon;
                        
                    elseif k1 == 1 && k2 == 3
                        epsilon = 10^-2;
                        f_tp_3 = friction_tp (initial_guess_alpha_g (j,1), initial_guess_alpha_w (j,1), rho_mixture_in1, initial_guess_V_m (j,1)+epsilon);
                        if j~=1
                            R1_3 = Residiual1 (initial_guess_P (j-1,1), initial_guess_P (j,1), rho_mixture_i1n1, rho_mixture_in, rho_mixture_in1, initial_guess_V_m (j+1,1)+epsilon , initial_V_m (j,1) + epsilon , initial_guess_V_m (j,1) + epsilon, f_tp_3 ) ;
                        else
                            R1_3 = Residiual1 (779100, initial_guess_P (j,1), rho_mixture_i1n1, rho_mixture_in, rho_mixture_in1, initial_guess_V_m (j+1,1) + epsilon , initial_V_m (j,1) + epsilon , initial_guess_V_m (j,1) + epsilon, f_tp_3 ) ;
                        end 
                        mat (k1, k2) = (R1_3 - R1_base) / epsilon;
                        
                    elseif k1 == 1 && k2 == 4 
                        epsilon = 10^2;
                        rho_water_in_4 = rho_P (14.7, (14.7*(initial_P (j,1)+epsilon)/101325), 10^-6, den_ref_water) * rho_T (20+273.15,  initial_T (j,1), 210*10^-6, den_ref_water) / den_ref_water;
                        rho_water_in1_4 = rho_P (14.7, (14.7*(initial_guess_P (j,1)+epsilon)/101325), 10^-6, den_ref_water) * rho_T (20+273.15, initial_guess_T (j,1), 210*10^-6, den_ref_water) / den_ref_water;
                        rho_water_i1n1_4 = rho_P (14.7, (14.7*(initial_guess_P (j+1,1)+epsilon)/101325), 10^-6, den_ref_water) * rho_T (20+273.15, initial_guess_T (j+1,1), 210*10^-6, den_ref_water) / den_ref_water;
                        rho_oil_in_4 = rho_P (14.7, (14.7*(initial_P (j,1)+epsilon)/101325), 6*10^-6, den_ref_oil) * rho_T (20+273.15,  initial_T (j,1), 950*10^-6, den_ref_oil) / den_ref_oil;
                        rho_oil_in1_4 = rho_P (14.7, (14.7*(initial_guess_P (j,1)+epsilon)/101325), 6*10^-6, den_ref_oil) * rho_T (20+273.15, initial_guess_T (j,1), 950*10^-6, den_ref_oil) / den_ref_oil;
                        rho_oil_i1n1_4 = rho_P (14.7, (14.7*(initial_guess_P (j+1,1)+epsilon)/101325), 6*10^-6, den_ref_oil) * rho_T (20+273.15, initial_guess_T (j+1,1), 950*10^-6, den_ref_oil) / den_ref_oil;                
                        rho_gas_in_4 = rho_P (14.7, (14.7*(initial_P (j,1)+epsilon)/101325), 100*10^-6, den_ref_gas) * rho_T (20+273.15,  initial_T (j,1), 3400*10^-6, den_ref_gas) / den_ref_gas;
                        rho_gas_in1_4 = rho_P (14.7, (14.7*(initial_guess_P (j,1)+epsilon)/101325), 100*10^-6, den_ref_gas) * rho_T (20+273.15, initial_guess_T (j,1), 3400*10^-6, den_ref_gas) / den_ref_gas;
                        rho_gas_i1n1_4 = rho_P (14.7, (14.7*(initial_guess_P (j+1,1)+epsilon)/101325), 100*10^-6, den_ref_gas) * rho_T (20+273.15, initial_guess_T (j+1,1), 3400*10^-6, den_ref_gas) / den_ref_gas;         
                        rho_mixture_in_4 = mixure_density (initial_alpha_gas (j,1), initial_alpha_water (j,1), rho_gas_in_4, rho_water_in_4,  rho_oil_in_4);
                        rho_mixture_in1_4 = mixure_density (initial_guess_alpha_g (j,1), initial_guess_alpha_w (j,1), rho_gas_in1_4, rho_water_in1_4,  rho_oil_in1_4);
                        rho_mixture_i1n1_4 = mixure_density (initial_guess_alpha_g (j+1,1), initial_guess_alpha_w (j+1,1),  rho_gas_i1n1_4, rho_water_i1n1_4,  rho_oil_i1n1_4);
                        f_tp_4 = friction_tp (initial_guess_alpha_g (j,1), initial_guess_alpha_w (j,1), rho_mixture_in1_4, initial_guess_V_m (j,1));
                        if j~=1
                            R1_4 = Residiual1 (initial_guess_P (j-1,1)+epsilon, initial_guess_P(j,1)+epsilon, rho_mixture_i1n1_4, rho_mixture_in_4, rho_mixture_in1_4, initial_guess_V_m (j+1,1) , initial_V_m (j,1), initial_guess_V_m (j,1), f_tp_4 ) ;
                        else
                            R1_4 = Residiual1 (779100+epsilon,initial_guess_P (j,1)+epsilon, rho_mixture_i1n1_4, rho_mixture_in_4, rho_mixture_in1_4, initial_guess_V_m (j+1,1) , initial_V_m (j,1), initial_guess_V_m (j,1), f_tp_4 ) ;
                        end                    
                        mat (k1, k2) = (R1_4 - R1_base) / epsilon;

                    elseif k1 == 1 && k2 == 5 
                        epsilon = 2;
                        rho_water_in_5 = rho_P (14.7, (14.7*initial_P (j,1)/101325), 10^-6, den_ref_water) * rho_T (20+273.15,  initial_T (j,1)+epsilon, 210*10^-6, den_ref_water) / den_ref_water;
                        rho_water_in1_5 = rho_P (14.7, (14.7*initial_guess_P (j,1)/101325), 10^-6, den_ref_water) * rho_T (20+273.15, initial_guess_T (j,1)+epsilon, 210*10^-6, den_ref_water) / den_ref_water;
                        rho_water_i1n1_5 = rho_P (14.7, (14.7*initial_guess_P (j+1,1)/101325), 10^-6, den_ref_water) * rho_T (20+273.15, initial_guess_T (j+1,1)+epsilon, 210*10^-6, den_ref_water) / den_ref_water;
                        rho_oil_in_5 = rho_P (14.7, (14.7*initial_P (j,1)/101325), 6*10^-6, den_ref_oil) * rho_T (20+273.15,  initial_T (j,1)+epsilon, 950*10^-6, den_ref_oil) / den_ref_oil;
                        rho_oil_in1_5 = rho_P (14.7, (14.7*initial_guess_P (j,1)/101325), 6*10^-6, den_ref_oil) * rho_T (20+273.15, initial_guess_T (j,1)+epsilon, 950*10^-6, den_ref_oil) / den_ref_oil;
                        rho_oil_i1n1_5 = rho_P (14.7, (14.7*initial_guess_P (j+1,1)/101325), 6*10^-6, den_ref_oil) * rho_T (20+273.15, initial_guess_T (j+1,1)+epsilon, 950*10^-6, den_ref_oil) / den_ref_oil;                    
                        rho_gas_in_5 = rho_P (14.7, (14.7*initial_P (j,1)/101325), 100*10^-6, den_ref_gas) * rho_T (20+273.15,  initial_T (j,1)+epsilon, 3400*10^-6, den_ref_gas) / den_ref_gas;
                        rho_gas_in1_5 = rho_P (14.7, (14.7*initial_guess_P (j,1)/101325), 100*10^-6, den_ref_gas) * rho_T (20+273.15, initial_guess_T (j,1)+epsilon, 3400*10^-6, den_ref_gas) / den_ref_gas;
                        rho_gas_i1n1_5 = rho_P (14.7, (14.7*initial_guess_P (j+1,1)/101325), 100*10^-6, den_ref_gas) * rho_T (20+273.15, initial_guess_T (j+1,1)+epsilon, 3400*10^-6, den_ref_gas) / den_ref_gas;         
                        rho_mixture_in_5 = mixure_density (initial_alpha_gas (j,1), initial_alpha_water (j,1), rho_gas_in_5, rho_water_in_5,  rho_oil_in_5);
                        rho_mixture_in1_5 = mixure_density (initial_guess_alpha_g (j,1), initial_guess_alpha_w (j,1), rho_gas_in1_5, rho_water_in1_5,  rho_oil_in1_5);
                        rho_mixture_i1n1_5 = mixure_density (initial_guess_alpha_g (j+1,1), initial_guess_alpha_w (j+1,1),  rho_gas_i1n1_5, rho_water_i1n1_5,  rho_oil_i1n1_5);
                        f_tp_5 = friction_tp (initial_guess_alpha_g (j,1), initial_guess_alpha_w (j,1), rho_mixture_in1_5, initial_guess_V_m (j,1));
                        if j~=1
                            R1_5 = Residiual1 (initial_guess_P (j-1,1), initial_guess_P (j,1), rho_mixture_i1n1_5, rho_mixture_in_5, rho_mixture_in1_5, initial_guess_V_m (j+1,1) , initial_V_m (j,1), initial_guess_V_m (j,1), f_tp ) ;
                        else
                            R1_5 = Residiual1 (779100, initial_guess_P (j,1), rho_mixture_i1n1_5, rho_mixture_in_5, rho_mixture_in1_5, initial_guess_V_m (j+1,1) , initial_V_m (j,1), initial_guess_V_m (j,1), f_tp ) ;
                        end  
                        mat (k1, k2) = (R1_5 - R1_base) / epsilon;

                   elseif k1 == 2 && k2 == 1
                        epsilon = 10^-2;
                        R2_1 = Residiual2 (rho_water_in, rho_water_in1, rho_water_i1n1, initial_alpha_water (j,1), initial_guess_alpha_w (j,1), initial_guess_alpha_w (j+1,1), q_water/Area_V, q_water);
                        mat (k1, k2) = (R2_1 - R2_base)/epsilon;
                        
                   elseif k1 == 2 && k2 == 2   
                        epsilon = 10^-2;
                        R2_2 = Residiual2 (rho_water_in, rho_water_in1, rho_water_i1n1, initial_alpha_water (j,1)+ epsilon, initial_guess_alpha_w (j,1)+ epsilon, initial_guess_alpha_w (j+1,1)+ epsilon, q_water/Area_V, q_water);
                        mat (k1, k2) = (R2_2 - R2_base)/epsilon;
                        
                   elseif k1 == 2 && k2 == 3
                        epsilon = 10^-2;
                        R2_3 = Residiual2 (rho_water_in, rho_water_in1, rho_water_i1n1, initial_alpha_water (j,1), initial_guess_alpha_w (j,1), initial_guess_alpha_w (j+1,1), q_water/Area_V, q_water);
                        mat (k1, k2) = (R2_3 - R2_base) / epsilon;

                   elseif k1 == 2 && k2 == 4
                        epsilon = 10^3;
                        rho_water_in_4 = rho_P (14.7, (14.7*(initial_P (j,1)+epsilon)/101325), 10^-6, den_ref_water) * rho_T (20+273.15,  initial_T (j,1), 210*10^-6, den_ref_water) / den_ref_water;
                        rho_water_in1_4 = rho_P (14.7, (14.7*(initial_guess_P (j,1)+epsilon)/101325), 10^-6, den_ref_water) * rho_T (20+273.15, initial_guess_T (j,1), 210*10^-6, den_ref_water) / den_ref_water;
                        rho_water_i1n1_4 = rho_P (14.7, (14.7*(initial_guess_P (j+1,1)+epsilon)/101325), 10^-6, den_ref_water) * rho_T (20+273.15, initial_guess_T (j+1,1), 210*10^-6, den_ref_water) / den_ref_water;
                        %=====
                        R2_4 = Residiual2 (rho_water_in_4, rho_water_in1_4, rho_water_i1n1_4, initial_alpha_water (j,1), initial_guess_alpha_w (j,1), initial_guess_alpha_w (j+1,1), q_water/Area_V, q_water);
                        mat (k1, k2) = (R2_4 - R2_base) / epsilon;

                   elseif k1 == 2 && k2 == 5
                        epsilon = 2;
                        rho_water_in_5 = rho_P (14.7, (14.7*initial_P (j,1)/101325), 10^-6, den_ref_water) * rho_T (20+273.15,  initial_T (j,1)+epsilon, 210*10^-6, den_ref_water) / den_ref_water;
                        rho_water_in1_5 = rho_P (14.7, (14.7*initial_guess_P (j,1)/101325), 10^-6, den_ref_water) * rho_T (20+273.15, initial_guess_T (j,1)+epsilon, 210*10^-6, den_ref_water) / den_ref_water;
                        rho_water_i1n1_5 = rho_P (14.7, (14.7*initial_guess_P (j+1,1)/101325), 10^-6, den_ref_water) * rho_T (20+273.15, initial_guess_T (j+1,1)+epsilon, 210*10^-6, den_ref_water) / den_ref_water;

                        R2_5 = Residiual2 (rho_water_in_5, rho_water_in1_5, rho_water_i1n1_5, initial_alpha_water (j,1), initial_guess_alpha_w (j,1), initial_guess_alpha_w (j+1,1), q_water/Area_V, q_water);
                        mat (k1, k2) = (R2_5 - R2_base) / epsilon;
                        
                   elseif k1 == 3 && k2 == 1
                        epsilon = 10^-2;
                        R3_1 = Residiual3 (rho_gas_in, rho_gas_in1,rho_gas_i1n1, initial_alpha_gas (j,1)+epsilon, initial_guess_alpha_g (j,1)+epsilon, initial_guess_alpha_g (j+1,1)+epsilon, ...
                        rho_oil_in, rho_oil_in1, rho_oil_i1n1,  initial_alpha_oil (j,1), initial_guess_alpha_o (j,1), initial_guess_alpha_o (j+1,1), 0, 1, q_oil, 0,...
                        gas_velocity_next_x_time, gas_velocity_next_time, initial_guess_V_m(j+1,1)-gas_velocity_next_x_time-q_water/Area_V, initial_guess_V_m(j,1)-gas_velocity_next_time-q_water/Area_V);
                        mat (k1, k2) = (R3_1 - R3_base) / epsilon;
                        
                   elseif k1 == 3 && k2 == 2
                        epsilon = 10^-2;
                        R3_2 = Residiual3 (rho_gas_in, rho_gas_in1,rho_gas_i1n1, initial_alpha_gas (j,1), initial_guess_alpha_g (j,1), initial_guess_alpha_g (j+1,1), ...
                        rho_oil_in, rho_oil_in1, rho_oil_i1n1,  initial_alpha_oil (j,1)-epsilon, initial_guess_alpha_o (j,1)-epsilon, initial_guess_alpha_o (j+1,1)-epsilon, 0, 1, q_oil, 0,...
                        gas_velocity_next_x_time, gas_velocity_next_time, initial_guess_V_m(j+1,1)-gas_velocity_next_x_time-q_water/Area_V, initial_guess_V_m(j,1)-gas_velocity_next_time-q_water/Area_V);
                        mat (k1, k2) = (R3_2 - R3_base) / epsilon;
                        
                   elseif k1 == 3 && k2 == 3
                        epsilon = 0.01;
                        gas_velocity_next_x_time_3 = gas_velocity (C0, V_d, initial_guess_V_m(j+1,1)+epsilon );
                        gas_velocity_next_time_3 = gas_velocity (C0, V_d, initial_guess_V_m(j,1)+epsilon);
                        R3_3 = Residiual3 (rho_gas_in, rho_gas_in1,rho_gas_i1n1, initial_alpha_gas (j,1), initial_guess_alpha_g (j,1), initial_guess_alpha_g (j+1,1), ...
                        rho_oil_in, rho_oil_in1, rho_oil_i1n1,  initial_alpha_oil (j,1), initial_guess_alpha_o (j,1), initial_guess_alpha_o (j+1,1), 0, 1, q_oil, 0,...
                        gas_velocity_next_x_time_3, gas_velocity_next_time_3, initial_guess_V_m(j+1,1)-gas_velocity_next_x_time_3-q_water/Area_V, initial_guess_V_m(j,1)-gas_velocity_next_time_3-q_water/Area_V);
                        mat (k1, k2) = (R3_3 - R3_base) / epsilon;
                        
                    elseif k1 == 3 && k2 == 4
                        epsilon = 10^3;
                        rho_water_in_4 = rho_P (14.7, (14.7*(initial_P (j,1)+epsilon)/101325), 10^-6, den_ref_water) * rho_T (20+273.15,  initial_T (j,1), 210*10^-6, den_ref_water) / den_ref_water;
                        rho_water_in1_4 = rho_P (14.7, (14.7*(initial_guess_P (j,1)+epsilon)/101325), 10^-6, den_ref_water) * rho_T (20+273.15, initial_guess_T (j,1), 210*10^-6, den_ref_water) / den_ref_water;
                        rho_water_i1n1_4 = rho_P (14.7, (14.7*(initial_guess_P (j+1,1)+epsilon)/101325), 10^-6, den_ref_water) * rho_T (20+273.15, initial_guess_T (j+1,1), 210*10^-6, den_ref_water) / den_ref_water;
                        rho_oil_in_4 = rho_P (14.7, (14.7*(initial_P (j,1)+epsilon)/101325), 6*10^-6, den_ref_oil) * rho_T (20+273.15,  initial_T (j,1), 950*10^-6, den_ref_oil) / den_ref_oil;
                        rho_oil_in1_4 = rho_P (14.7, (14.7*(initial_guess_P (j,1)+epsilon)/101325), 6*10^-6, den_ref_oil) * rho_T (20+273.15, initial_guess_T (j,1), 950*10^-6, den_ref_oil) / den_ref_oil;
                        rho_oil_i1n1_4 = rho_P (14.7, (14.7*(initial_guess_P (j+1,1)+epsilon)/101325), 6*10^-6, den_ref_oil) * rho_T (20+273.15, initial_guess_T (j+1,1), 950*10^-6, den_ref_oil) / den_ref_oil;      
                        rho_gas_in_4 = rho_P (14.7, (14.7*(initial_P (j,1)+epsilon)/101325), 100*10^-6, den_ref_gas) * rho_T (20+273.15,  initial_T (j,1), 3400*10^-6, den_ref_gas) / den_ref_gas;
                        rho_gas_in1_4 = rho_P (14.7, (14.7*(initial_guess_P (j,1)+epsilon)/101325), 100*10^-6, den_ref_gas) * rho_T (20+273.15, initial_guess_T (j,1), 3400*10^-6, den_ref_gas) / den_ref_gas;
                        rho_gas_i1n1_4 = rho_P (14.7, (14.7*(initial_guess_P (j+1,1)+epsilon)/101325), 100*10^-6, den_ref_gas) * rho_T (20+273.15, initial_guess_T (j+1,1), 3400*10^-6, den_ref_gas) / den_ref_gas;         
                        R3_4 = Residiual3 (rho_gas_in_4, rho_gas_in1_4,rho_gas_i1n1_4, initial_alpha_gas (j,1), initial_guess_alpha_g (j,1), initial_guess_alpha_g (j+1,1), ...
                        rho_oil_in_4, rho_oil_in1_4, rho_oil_i1n1_4,  initial_alpha_oil (j,1), initial_guess_alpha_o (j,1), initial_guess_alpha_o (j+1,1), 0, 1, q_oil, 0,...
                        gas_velocity_next_x_time, gas_velocity_next_time, initial_guess_V_m(j+1,1)-gas_velocity_next_x_time-q_water/Area_V, initial_guess_V_m(j,1)-gas_velocity_next_time-q_water/Area_V);
                        mat (k1, k2) = (R3_4 - R3_base) / epsilon;
                        
                    elseif k1 == 3 && k2 == 5
                        epsilon = 2;
                        rho_water_in_5 = rho_P (14.7, (14.7*initial_P (j,1)/101325), 10^-6, den_ref_water) * rho_T (20+273.15,  initial_T (j,1)+epsilon, 210*10^-6, den_ref_water) / den_ref_water;
                        rho_water_in1_5 = rho_P (14.7, (14.7*initial_guess_P (j,1)/101325), 10^-6, den_ref_water) * rho_T (20+273.15, initial_guess_T (j,1)+epsilon, 210*10^-6, den_ref_water) / den_ref_water;
                        rho_water_i1n1_5 = rho_P (14.7, (14.7*initial_guess_P (j+1,1)/101325), 10^-6, den_ref_water) * rho_T (20+273.15, initial_guess_T (j+1,1)+epsilon, 210*10^-6, den_ref_water) / den_ref_water;
                        rho_oil_in_5 = rho_P (14.7, (14.7*initial_P (j,1)/101325), 6*10^-6, den_ref_oil) * rho_T (20+273.15,  initial_T (j,1)+epsilon, 950*10^-6, den_ref_oil) / den_ref_oil;
                        rho_oil_in1_5 = rho_P (14.7, (14.7*initial_guess_P (j,1)/101325), 6*10^-6, den_ref_oil) * rho_T (20+273.15, initial_guess_T (j,1)+epsilon, 950*10^-6, den_ref_oil) / den_ref_oil;
                        rho_oil_i1n1_5 = rho_P (14.7, (14.7*initial_guess_P (j+1,1)/101325), 6*10^-6, den_ref_oil) * rho_T (20+273.15, initial_guess_T (j+1,1)+epsilon, 950*10^-6, den_ref_oil) / den_ref_oil;      
                        rho_gas_in_5 = rho_P (14.7, (14.7*initial_P (j,1)/101325), 100*10^-6, den_ref_gas) * rho_T (20+273.15,  initial_T (j,1)+epsilon, 3400*10^-6, den_ref_gas) / den_ref_gas;
                        rho_gas_in1_5 = rho_P (14.7, (14.7*initial_guess_P (j,1)/101325), 100*10^-6, den_ref_gas) * rho_T (20+273.15, initial_guess_T (j,1)+epsilon, 3400*10^-6, den_ref_gas) / den_ref_gas;
                        rho_gas_i1n1_5 = rho_P (14.7, (14.7*initial_guess_P (j+1,1)/101325), 100*10^-6, den_ref_gas) * rho_T (20+273.15, initial_guess_T (j+1,1)+epsilon, 3400*10^-6, den_ref_gas) / den_ref_gas;         
                        R3_5 = Residiual3 (rho_gas_in_5, rho_gas_in1_5,rho_gas_i1n1_5, initial_alpha_gas (j,1), initial_guess_alpha_g (j,1), initial_guess_alpha_g (j+1,1), ...
                        rho_oil_in_5, rho_oil_in1_5, rho_oil_i1n1_5,  initial_alpha_oil (j,1), initial_guess_alpha_o (j,1), initial_guess_alpha_o (j+1,1), 0, 1, q_oil, 0,...
                        gas_velocity_next_x_time, gas_velocity_next_time, initial_guess_V_m(j+1,1)-gas_velocity_next_x_time-q_water/Area_V, initial_guess_V_m(j,1)-gas_velocity_next_time-q_water/Area_V);
                        mat (k1, k2) = (R3_5 - R3_base) / epsilon;
                        
                   elseif k1 == 4 && k2 == 1   
                       epsilon = 10^-2;
                       R4_1 = Residiual4 (rho_gas_in, rho_gas_in1, rho_gas_i1n1, initial_alpha_gas (j,1)+epsilon, initial_guess_alpha_g (j,1)+epsilon, initial_guess_alpha_g (j+1,1)+epsilon, ...
                        rho_oil_in, rho_oil_in1, rho_oil_i1n1, initial_alpha_oil (j,1), initial_guess_alpha_o (j,1), initial_guess_alpha_o (j+1,1),...
                        rho_water_in, rho_water_in1, rho_water_i1n1, initial_alpha_water (j,1), initial_guess_alpha_w (j,1), initial_guess_alpha_w (j+1,1), 1, 0, 0, q_gas, 0, 0,...
                        gas_velocity_next_x_time, gas_velocity_next_time, initial_guess_V_m(j+1,1)-gas_velocity_next_x_time-q_water/Area_V, initial_guess_V_m(j,1)-gas_velocity_next_time-q_water/Area_V, q_water/Area_V, q_water/Area_V) ;
                        mat (k1, k2) = (R4_1 - R4_base) / epsilon;
                        
                   elseif k1 == 4 && k2 == 2   
                       epsilon = 10^-2;
                       R4_2 = Residiual4 (rho_gas_in, rho_gas_in1, rho_gas_i1n1, initial_alpha_gas (j,1), initial_guess_alpha_g (j,1), initial_guess_alpha_g (j+1,1), ...
                        rho_oil_in, rho_oil_in1, rho_oil_i1n1, initial_alpha_oil (j,1), initial_guess_alpha_o (j,1), initial_guess_alpha_o (j+1,1),...
                        rho_water_in, rho_water_in1, rho_water_i1n1, initial_alpha_water (j,1)+epsilon, initial_guess_alpha_w (j,1)+epsilon, initial_guess_alpha_w (j+1,1)+epsilon, 1, 0, 0, q_gas, 0, 0,...
                        gas_velocity_next_x_time, gas_velocity_next_time, initial_guess_V_m(j+1,1)-gas_velocity_next_x_time-q_water/Area_V, initial_guess_V_m(j,1)-gas_velocity_next_time-q_water/Area_V, q_water/Area_V, q_water/Area_V) ;
                        mat (k1, k2) = (R4_2 - R4_base) / epsilon;
                        
                    elseif k1 == 4 && k2 == 3
                        epsilon = 0.01;
                        gas_velocity_next_x_time_3 = gas_velocity (C0, V_d, initial_guess_V_m(j+1,1)+epsilon );
                        gas_velocity_next_time_3 = gas_velocity (C0, V_d, initial_guess_V_m(j,1)+epsilon);
                         R4_3 = Residiual4 (rho_gas_in, rho_gas_in1, rho_gas_i1n1, initial_alpha_gas (j,1), initial_guess_alpha_g (j,1), initial_guess_alpha_g (j+1,1), ...
                        rho_oil_in, rho_oil_in1, rho_oil_i1n1, initial_alpha_oil (j,1), initial_guess_alpha_o (j,1), initial_guess_alpha_o (j+1,1),...
                        rho_water_in, rho_water_in1, rho_water_i1n1, initial_alpha_water (j,1)+epsilon, initial_guess_alpha_w (j,1)+epsilon, initial_guess_alpha_w (j+1,1)+epsilon, 1, 0, 0, q_gas, 0, 0,...
                        gas_velocity_next_x_time_3, gas_velocity_next_time_3, initial_guess_V_m(j+1,1)-gas_velocity_next_x_time_3-q_water/Area_V, initial_guess_V_m(j,1)-gas_velocity_next_time_3-q_water/Area_V, q_water/Area_V, q_water/Area_V) ;
                        mat (k1, k2) = (R4_3 - R4_base) / epsilon;    
                        
                    elseif k1 == 4 && k2 == 4
                        epsilon = 10^3;
                        rho_water_in_4 = rho_P (14.7, (14.7*(initial_P (j,1)+epsilon)/101325), 10^-6, den_ref_water) * rho_T (20+273.15,  initial_T (j,1), 210*10^-6, den_ref_water) / den_ref_water;
                        rho_water_in1_4 = rho_P (14.7, (14.7*(initial_guess_P (j,1)+epsilon)/101325), 10^-6, den_ref_water) * rho_T (20+273.15, initial_guess_T (j,1), 210*10^-6, den_ref_water) / den_ref_water;
                        rho_water_i1n1_4 = rho_P (14.7, (14.7*(initial_guess_P (j+1,1)+epsilon)/101325), 10^-6, den_ref_water) * rho_T (20+273.15, initial_guess_T (j+1,1), 210*10^-6, den_ref_water) / den_ref_water;
                        rho_oil_in_4 = rho_P (14.7, (14.7*(initial_P (j,1)+epsilon)/101325), 6*10^-6, den_ref_oil) * rho_T (20+273.15,  initial_T (j,1), 950*10^-6, den_ref_oil) / den_ref_oil;
                        rho_oil_in1_4 = rho_P (14.7, (14.7*(initial_guess_P (j,1)+epsilon)/101325), 6*10^-6, den_ref_oil) * rho_T (20+273.15, initial_guess_T (j,1), 950*10^-6, den_ref_oil) / den_ref_oil;
                        rho_oil_i1n1_4 = rho_P (14.7, (14.7*(initial_guess_P (j+1,1)+epsilon)/101325), 6*10^-6, den_ref_oil) * rho_T (20+273.15, initial_guess_T (j+1,1), 950*10^-6, den_ref_oil) / den_ref_oil;      
                        rho_gas_in_4 = rho_P (14.7, (14.7*(initial_P (j,1)+epsilon)/101325), 100*10^-6, den_ref_gas) * rho_T (20+273.15,  initial_T (j,1), 3400*10^-6, den_ref_gas) / den_ref_gas;
                        rho_gas_in1_4 = rho_P (14.7, (14.7*(initial_guess_P (j,1)+epsilon)/101325), 100*10^-6, den_ref_gas) * rho_T (20+273.15, initial_guess_T (j,1), 3400*10^-6, den_ref_gas) / den_ref_gas;
                        rho_gas_i1n1_4 = rho_P (14.7, (14.7*(initial_guess_P (j+1,1)+epsilon)/101325), 100*10^-6, den_ref_gas) * rho_T (20+273.15, initial_guess_T (j+1,1), 3400*10^-6, den_ref_gas) / den_ref_gas;         
                        R4_4 = Residiual4 (rho_gas_in_4, rho_gas_in1_4, rho_gas_i1n1_4, initial_alpha_gas (j,1), initial_guess_alpha_g (j,1), initial_guess_alpha_g (j+1,1), ...
                        rho_oil_in_4, rho_oil_in1_4, rho_oil_i1n1_4, initial_alpha_oil (j,1), initial_guess_alpha_o (j,1), initial_guess_alpha_o (j+1,1),...
                        rho_water_in_4, rho_water_in1_4, rho_water_i1n1_4, initial_alpha_water (j,1), initial_guess_alpha_w (j,1), initial_guess_alpha_w (j+1,1), 1, 0, 0, q_gas, 0, 0,...
                        gas_velocity_next_x_time, gas_velocity_next_time, initial_guess_V_m(j+1,1)-gas_velocity_next_x_time-q_water/Area_V, initial_guess_V_m(j,1)-gas_velocity_next_time-q_water/Area_V, q_water/Area_V, q_water/Area_V) ;
                        mat (k1, k2) = (R4_4 - R4_base) / epsilon; 
                        
                     elseif k1 == 4 && k2 == 5
                        epsilon = 2;
                        rho_water_in_5 = rho_P (14.7, (14.7*initial_P (j,1)/101325), 10^-6, den_ref_water) * rho_T (20+273.15,  initial_T (j,1)+epsilon, 210*10^-6, den_ref_water) / den_ref_water;
                        rho_water_in1_5 = rho_P (14.7, (14.7*initial_guess_P (j,1)/101325), 10^-6, den_ref_water) * rho_T (20+273.15, initial_guess_T (j,1)+epsilon, 210*10^-6, den_ref_water) / den_ref_water;
                        rho_water_i1n1_5 = rho_P (14.7, (14.7*initial_guess_P (j+1,1)/101325), 10^-6, den_ref_water) * rho_T (20+273.15, initial_guess_T (j+1,1)+epsilon, 210*10^-6, den_ref_water) / den_ref_water;
                        rho_oil_in_5 = rho_P (14.7, (14.7*initial_P (j,1)/101325), 6*10^-6, den_ref_oil) * rho_T (20+273.15,  initial_T (j,1)+epsilon, 950*10^-6, den_ref_oil) / den_ref_oil;
                        rho_oil_in1_5 = rho_P (14.7, (14.7*initial_guess_P (j,1)/101325), 6*10^-6, den_ref_oil) * rho_T (20+273.15, initial_guess_T (j,1)+epsilon, 950*10^-6, den_ref_oil) / den_ref_oil;
                        rho_oil_i1n1_5 = rho_P (14.7, (14.7*initial_guess_P (j+1,1)/101325), 6*10^-6, den_ref_oil) * rho_T (20+273.15, initial_guess_T (j+1,1)+epsilon, 950*10^-6, den_ref_oil) / den_ref_oil;                   
                        rho_gas_in_5 = rho_P (14.7, (14.7*initial_P (j,1)/101325), 100*10^-6, den_ref_gas) * rho_T (20+273.15,  initial_T (j,1)+epsilon, 3400*10^-6, den_ref_gas) / den_ref_gas;
                        rho_gas_in1_5 = rho_P (14.7, (14.7*initial_guess_P (j,1)/101325), 100*10^-6, den_ref_gas) * rho_T (20+273.15, initial_guess_T (j,1)+epsilon, 3400*10^-6, den_ref_gas) / den_ref_gas;
                        rho_gas_i1n1_5 = rho_P (14.7, (14.7*initial_guess_P (j+1,1)/101325), 100*10^-6, den_ref_gas) * rho_T (20+273.15, initial_guess_T (j+1,1)+epsilon, 3400*10^-6, den_ref_gas) / den_ref_gas;         
                        R4_5 = Residiual4 (rho_gas_in_5, rho_gas_in1_5, rho_gas_i1n1_5, initial_alpha_gas (j,1), initial_guess_alpha_g (j,1), initial_guess_alpha_g (j+1,1), ...
                        rho_oil_in_5, rho_oil_in1_5, rho_oil_i1n1_5, initial_alpha_oil (j,1), initial_guess_alpha_o (j,1), initial_guess_alpha_o (j+1,1),...
                        rho_water_in_5, rho_water_in1_5, rho_water_i1n1_5, initial_alpha_water (j,1), initial_guess_alpha_w (j,1), initial_guess_alpha_w (j+1,1), 1, 0, 0, q_gas, 0, 0,...
                        gas_velocity_next_x_time, gas_velocity_next_time, initial_guess_V_m(j+1,1)-gas_velocity_next_x_time-q_water/Area_V, initial_guess_V_m(j,1)-gas_velocity_next_time-q_water/Area_V, q_water/Area_V, q_water/Area_V) ;
                        mat (k1, k2) = (R4_5 - R4_base) / epsilon;
                        
                   elseif k1 == 5 && k2 == 1   
                        epsilon = 10^-2;
                        R5_1 = Residiual5 (rho_gas_in, rho_gas_in1, rho_gas_i1n1, initial_alpha_gas (j,1)+epsilon, initial_guess_alpha_g (j,1)+epsilon, initial_guess_alpha_g (j+1,1)+epsilon, ...
                        rho_oil_in, rho_oil_in1, rho_oil_i1n1, initial_alpha_oil (j,1), initial_guess_alpha_o (j,1), initial_guess_alpha_o (j+1,1),...
                        rho_water_in, rho_water_in1, rho_water_i1n1, initial_alpha_water (j,1), initial_guess_alpha_w (j,1), initial_guess_alpha_w (j+1,1),...
                        gas_velocity_next_x_time, gas_velocity_next_time, gas_velocity_present, initial_guess_V_m(j+1,1)-gas_velocity_next_x_time-q_water/Area, initial_guess_V_m(j,1)-gas_velocity_next_time-q_water/Area, initial_guess_V_m(j,1)-gas_velocity_present-q_water/Area,q_water/Area, q_water/Area, q_water/Area, ...
                        Up (Cp_G, initial_T (j,1), 20+273.15), Up (Cp_G, initial_guess_T (j,1), 20+273.15), Up (Cp_G, initial_guess_T (j+1,1), 20+273.15),...
                        Up (Cp_O, initial_T (j,1), 20+273.15), Up (Cp_O, initial_guess_T (j,1), 20+273.15), Up (Cp_O, initial_guess_T (j+1,1), 20+273.15),...
                        Up (Cp_W, initial_T (j,1), 20+273.15), Up (Cp_W, initial_guess_T (j,1), 20+273.15), Up (Cp_W, initial_guess_T (j+1,1), 20+273.15),...
                        initial_guess_P (j,1), initial_guess_P (j+1,1), Q_loss (T_top, T_bottom, U_to_top, U_to_end, initial_guess_T (j,1), j));
                        mat (k1, k2) = (R5_1 - R5_base) / epsilon;
                        
                    elseif k1 == 5 && k2 == 2   
                        epsilon = 10^-2;
                        R5_2 = Residiual5 (rho_gas_in, rho_gas_in1, rho_gas_i1n1, initial_alpha_gas (j,1), initial_guess_alpha_g (j,1), initial_guess_alpha_g (j+1,1), ...
                        rho_oil_in, rho_oil_in1, rho_oil_i1n1, initial_alpha_oil (j,1), initial_guess_alpha_o (j,1), initial_guess_alpha_o (j+1,1),...
                        rho_water_in, rho_water_in1, rho_water_i1n1, initial_alpha_water (j,1)+epsilon, initial_guess_alpha_w (j,1)+epsilon, initial_guess_alpha_w (j+1,1)+epsilon,...
                        gas_velocity_next_x_time, gas_velocity_next_time, gas_velocity_present, initial_guess_V_m(j+1,1)-gas_velocity_next_x_time-q_water/Area, initial_guess_V_m(j,1)-gas_velocity_next_time-q_water/Area, initial_guess_V_m(j,1)-gas_velocity_present-q_water/Area,q_water/Area, q_water/Area, q_water/Area, ...
                        Up (Cp_G, initial_T (j,1), 20+273.15), Up (Cp_G, initial_guess_T (j,1), 20+273.15), Up (Cp_G, initial_guess_T (j+1,1), 20+273.15),...
                        Up (Cp_O, initial_T (j,1), 20+273.15), Up (Cp_O, initial_guess_T (j,1), 20+273.15), Up (Cp_O, initial_guess_T (j+1,1), 20+273.15),...
                        Up (Cp_W, initial_T (j,1), 20+273.15), Up (Cp_W, initial_guess_T (j,1), 20+273.15), Up (Cp_W, initial_guess_T (j+1,1), 20+273.15),...
                        initial_guess_P (j,1), initial_guess_P (j+1,1), Q_loss (T_top, T_bottom, U_to_top, U_to_end, initial_guess_T (j,1), j));
                        mat (k1, k2) = (R5_2 - R5_base) / epsilon ;   
                        
                   elseif k1 == 5 && k2 == 3   
                        epsilon = 10^-2;                  
                        gas_velocity_present_3 = gas_velocity (C0, V_d, initial_V_m(j,1));
                        gas_velocity_next_time_3 = gas_velocity (C0, V_d, initial_guess_V_m(j,1));
                        gas_velocity_next_x_time_3 = gas_velocity (C0, V_d, initial_guess_V_m(j+1,1));
                        R5_3 = Residiual5 (rho_gas_in, rho_gas_in1, rho_gas_i1n1, initial_alpha_gas (j,1), initial_guess_alpha_g (j,1), initial_guess_alpha_g (j+1,1), ...
                        rho_oil_in, rho_oil_in1, rho_oil_i1n1, initial_alpha_oil (j,1), initial_guess_alpha_o (j,1), initial_guess_alpha_o (j+1,1),...
                        rho_water_in, rho_water_in1, rho_water_i1n1, initial_alpha_water (j,1)+epsilon, initial_guess_alpha_w (j,1)+epsilon, initial_guess_alpha_w (j+1,1)+epsilon,...
                        gas_velocity_next_x_time_3, gas_velocity_next_time_3, gas_velocity_present_3, initial_guess_V_m(j+1,1)-gas_velocity_next_x_time_3-q_water/Area, initial_guess_V_m(j,1)-gas_velocity_next_time_3-q_water/Area, initial_guess_V_m(j,1)-gas_velocity_present_3-q_water/Area,q_water/Area, q_water/Area, q_water/Area, ...
                        Up (Cp_G, initial_T (j,1), 20+273.15), Up (Cp_G, initial_guess_T (j,1), 20+273.15), Up (Cp_G, initial_guess_T (j+1,1), 20+273.15),...
                        Up (Cp_O, initial_T (j,1), 20+273.15), Up (Cp_O, initial_guess_T (j,1), 20+273.15), Up (Cp_O, initial_guess_T (j+1,1), 20+273.15),...
                        Up (Cp_W, initial_T (j,1), 20+273.15), Up (Cp_W, initial_guess_T (j,1), 20+273.15), Up (Cp_W, initial_guess_T (j+1,1), 20+273.15),...
                        initial_guess_P (j,1), initial_guess_P (j+1,1), Q_loss (T_top, T_bottom, U_to_top, U_to_end, initial_guess_T (j,1), j));
                        mat (k1, k2) = (R5_3 - R5_base) / epsilon ;  
                        
                  elseif k1 == 5 && k2 == 4 
                        epsilon = 10^3;
                        rho_water_in_4 = rho_P (14.7, (14.7*(initial_P (j,1)+epsilon)/101325), 10^-6, den_ref_water) * rho_T (20+273.15,  initial_T (j,1), 210*10^-6, den_ref_water) / den_ref_water;
                        rho_water_in1_4 = rho_P (14.7, (14.7*(initial_guess_P (j,1)+epsilon)/101325), 10^-6, den_ref_water) * rho_T (20+273.15, initial_guess_T (j,1), 210*10^-6, den_ref_water) / den_ref_water;
                        rho_water_i1n1_4 = rho_P (14.7, (14.7*(initial_guess_P (j+1,1)+epsilon)/101325), 10^-6, den_ref_water) * rho_T (20+273.15, initial_guess_T (j+1,1), 210*10^-6, den_ref_water) / den_ref_water;
                        rho_oil_in_4 = rho_P (14.7, (14.7*(initial_P (j,1)+epsilon)/101325), 6*10^-6, den_ref_oil) * rho_T (20+273.15,  initial_T (j,1), 950*10^-6, den_ref_oil) / den_ref_oil;
                        rho_oil_in1_4 = rho_P (14.7, (14.7*(initial_guess_P (j,1)+epsilon)/101325), 6*10^-6, den_ref_oil) * rho_T (20+273.15, initial_guess_T (j,1), 950*10^-6, den_ref_oil) / den_ref_oil;
                        rho_oil_i1n1_4 = rho_P (14.7, (14.7*(initial_guess_P (j+1,1)+epsilon)/101325), 6*10^-6, den_ref_oil) * rho_T (20+273.15, initial_guess_T (j+1,1), 950*10^-6, den_ref_oil) / den_ref_oil;                   
                        rho_gas_in_4 = rho_P (14.7, (14.7*(initial_P (j,1)+epsilon)/101325), 100*10^-6, den_ref_gas) * rho_T (20+273.15,  initial_T (j,1), 3400*10^-6, den_ref_gas) / den_ref_gas;
                        rho_gas_in1_4 = rho_P (14.7, (14.7*(initial_guess_P (j,1)+epsilon)/101325), 100*10^-6, den_ref_gas) * rho_T (20+273.15, initial_guess_T (j,1), 3400*10^-6, den_ref_gas) / den_ref_gas;
                        rho_gas_i1n1_4 = rho_P (14.7, (14.7*(initial_guess_P (j+1,1)+epsilon)/101325), 100*10^-6, den_ref_gas) * rho_T (20+273.15, initial_guess_T (j+1,1), 3400*10^-6, den_ref_gas) / den_ref_gas;         
                        
                        R5_4 = Residiual5 (rho_gas_in_4, rho_gas_in1_4, rho_gas_i1n1_4, initial_alpha_gas (j,1), initial_guess_alpha_g (j,1), initial_guess_alpha_g (j+1,1), ...
                        rho_oil_in_4, rho_oil_in1_4, rho_oil_i1n1_4, initial_alpha_oil (j,1), initial_guess_alpha_o (j,1), initial_guess_alpha_o (j+1,1),...
                        rho_water_in_4, rho_water_in1_4, rho_water_i1n1_4, initial_alpha_water (j,1), initial_guess_alpha_w (j,1), initial_guess_alpha_w (j+1,1),...
                        gas_velocity_next_x_time, gas_velocity_next_time, gas_velocity_present, initial_guess_V_m(j+1,1)-gas_velocity_next_x_time-q_water/Area, initial_guess_V_m(j,1)-gas_velocity_next_time-q_water/Area, initial_guess_V_m(j,1)-gas_velocity_present-q_water/Area,q_water/Area, q_water/Area, q_water/Area, ...
                        Up (Cp_G, initial_T (j,1), 20+273.15), Up (Cp_G, initial_guess_T (j,1), 20+273.15), Up (Cp_G, initial_guess_T (j+1,1), 20+273.15),...
                        Up (Cp_O, initial_T (j,1), 20+273.15), Up (Cp_O, initial_guess_T (j,1), 20+273.15), Up (Cp_O, initial_guess_T (j+1,1), 20+273.15),...
                        Up (Cp_W, initial_T (j,1), 20+273.15), Up (Cp_W, initial_guess_T (j,1), 20+273.15), Up (Cp_W, initial_guess_T (j+1,1), 20+273.15),...
                        initial_guess_P (j,1)+epsilon, initial_guess_P (j+1,1)+epsilon, Q_loss (T_top, T_bottom, U_to_top, U_to_end, initial_guess_T (j,1), j));
                        mat (k1, k2) = (R5_4 - R5_base) / epsilon ;
                        
                    elseif k1 == 5 && k2 == 5 
                        epsilon = 2;
                        rho_water_in_5 = rho_P (14.7, (14.7*initial_P (j,1)/101325), 10^-6, den_ref_water) * rho_T (20+273.15,  initial_T (j,1)+epsilon, 210*10^-6, den_ref_water) / den_ref_water;
                        rho_water_in1_5 = rho_P (14.7, (14.7*initial_guess_P (j,1)/101325), 10^-6, den_ref_water) * rho_T (20+273.15, initial_guess_T (j,1)+epsilon, 210*10^-6, den_ref_water) / den_ref_water;
                        rho_water_i1n1_5 = rho_P (14.7, (14.7*initial_guess_P (j+1,1)/101325), 10^-6, den_ref_water) * rho_T (20+273.15, initial_guess_T (j+1,1)+epsilon, 210*10^-6, den_ref_water) / den_ref_water;
                        rho_oil_in_5 = rho_P (14.7, (14.7*initial_P (j,1)/101325), 6*10^-6, den_ref_oil) * rho_T (20+273.15,  initial_T (j,1)+epsilon, 950*10^-6, den_ref_oil) / den_ref_oil;
                        rho_oil_in1_5 = rho_P (14.7, (14.7*initial_guess_P (j,1)/101325), 6*10^-6, den_ref_oil) * rho_T (20+273.15, initial_guess_T (j,1)+epsilon, 950*10^-6, den_ref_oil) / den_ref_oil;
                        rho_oil_i1n1_5 = rho_P (14.7, (14.7*initial_guess_P (j+1,1)/101325), 6*10^-6, den_ref_oil) * rho_T (20+273.15, initial_guess_T (j+1,1)+epsilon, 950*10^-6, den_ref_oil) / den_ref_oil;               
                        rho_gas_in_5 = rho_P (14.7, (14.7*initial_P (j,1)/101325), 100*10^-6, den_ref_gas) * rho_T (20+273.15,  initial_T (j,1)+epsilon, 3400*10^-6, den_ref_gas) / den_ref_gas;
                        rho_gas_in1_5 = rho_P (14.7, (14.7*initial_guess_P (j,1)/101325), 100*10^-6, den_ref_gas) * rho_T (20+273.15, initial_guess_T (j,1)+epsilon, 3400*10^-6, den_ref_gas) / den_ref_gas;
                        rho_gas_i1n1_5 = rho_P (14.7, (14.7*initial_guess_P (j+1,1)/101325), 100*10^-6, den_ref_gas) * rho_T (20+273.15, initial_guess_T (j+1,1)+epsilon, 3400*10^-6, den_ref_gas) / den_ref_gas;         
                        
                        R5_5 = Residiual5 (rho_gas_in_5, rho_gas_in1_5, rho_gas_i1n1_5, initial_alpha_gas (j,1), initial_guess_alpha_g (j,1), initial_guess_alpha_g (j+1,1), ...
                        rho_oil_in_5, rho_oil_in1_5, rho_oil_i1n1_5, initial_alpha_oil (j,1), initial_guess_alpha_o (j,1), initial_guess_alpha_o (j+1,1),...
                        rho_water_in_5, rho_water_in1_5, rho_water_i1n1_5, initial_alpha_water (j,1), initial_guess_alpha_w (j,1), initial_guess_alpha_w (j+1,1),...
                        gas_velocity_next_x_time, gas_velocity_next_time, gas_velocity_present, initial_guess_V_m(j+1,1)-gas_velocity_next_x_time-q_water/Area, initial_guess_V_m(j,1)-gas_velocity_next_time-q_water/Area, initial_guess_V_m(j,1)-gas_velocity_present-q_water/Area,q_water/Area, q_water/Area, q_water/Area, ...
                        Up (Cp_G, initial_T (j,1)+epsilon, 20+273.15), Up (Cp_G, initial_guess_T (j,1)+epsilon, 20+273.15), Up (Cp_G, initial_guess_T (j+1,1)+epsilon, 20+273.15),...
                        Up (Cp_O, initial_T (j,1)+epsilon, 20+273.15), Up (Cp_O, initial_guess_T (j,1)+epsilon, 20+273.15), Up (Cp_O, initial_guess_T (j+1,1)+epsilon, 20+273.15),...
                        Up (Cp_W, initial_T (j,1)+epsilon, 20+273.15), Up (Cp_W, initial_guess_T (j,1)+epsilon, 20+273.15), Up (Cp_W, initial_guess_T (j+1,1)+epsilon, 20+273.15),...
                        initial_guess_P (j,1), initial_guess_P (j+1,1), Q_loss (T_top, T_bottom, U_to_top, U_to_end, initial_guess_T (j,1)+epsilon, j));
                        mat (k1, k2) = (R5_5 - R5_base) / epsilon;
                    end
                end
            end

        cell1{j} = mat;
        clear mat
        end
        for k = 1: number_grids 
            if k == 1
                t = 5;
                Jacobian (1:t,1:5) = cell1{1};
                Jacobian (1:t,6:10) = cell1{2};
                answer (1:t,1) = - cell2{1};
            elseif k == 2
                Jacobian (t+1:t+5,1:5) = cell1{1};
                Jacobian (t+1:t+5,6:10) = cell1{2};
                Jacobian (t+1:t+5,11:15) = cell1{3};
                answer (t+1:t+5,1) = -cell2{2};
                t = t+5;
            elseif k > 2 && k<= number_grids-1
                Jacobian (t+1:t+5,(k-2)*5+1:((k-2)*5+5)) = cell1{k-1};
                a = (k-2)*5+5;
                Jacobian (t+1:t+5,a+1:a+5) = cell1{k};
                b = a+5;
                Jacobian (t+1:t+5,b+1:b+5) = cell1{k+1};  
                answer (t+1:t+5,1) = -cell2{k};
                t = t+5;
            elseif k== number_grids
                Jacobian (t+1:t+5,(k-2)*5+1:((k-2)*5+5)) = cell1{k-1};
                a = (k-2)*5+5;
                Jacobian (t+1:t+5,a+1:a+5) = cell1{k};
                answer (t+1:t+5,1) = -cell2{k};             
            end
        end
        initial_guess = zeros (length(answer),1);
        update = SOR_method (Jacobian,-answer,1.5, initial_guess, 10^-2); 
        alpha_g2 = initial_guess_alpha_g (1:end-1) + update(1:5:end);
        alpha_w2 = initial_guess_alpha_w (1:end-1)+ update(2:5:end,1);
        V_m2 = initial_guess_V_m (1:end-1)+ update(3:5:end,1);
        pressure2 = initial_guess_P (1:end-1)+ update(4:5:end,1);
        temp2 = initial_guess_T (1:end-1) + update(5:5:end,1);
        Flag = true;
        
        if update(1:5:end) > 10^-4
            Flag = false;
        end
        
        if update(2:5:end) > 10^-4
            Flag = false;
        end
        
        if update(3:5:end) > 10^-4
            Flag = false;
        end
        
        if update(4:5:end) > 10
            Flag = false;
        end
        
        if update(4:5:end) > 2
            Flag = false;
        end
        
        if Flag == true % next timestep
            initial_T = temp2;
            initial_P = pressure2; 
            initial_alpha_gas = alpha_g2;
            initial_alpha_water  = alpha_w2;
            initial_V_m = V_m2;
            initial_alpha_oil = 1-initial_alpha_gas-initial_alpha_water;
            break
            
        else % more iteration
            initial_guess_alpha_g = alpha_g2;
            initial_guess_alpha_w = alpha_w2;
            initial_guess_alpha_o = 1-  initial_guess_alpha_g(:,1)- initial_guess_alpha_w(:,1);
            initial_guess_V_m = V_m2;
             initial_V_m (1,1) = 0;
            initial_guess_P = pressure2;
            initial_guess_T =  temp2;
        end
    end            
end

x = [0: number_grids]';
y=temp2-273.15;
y=[31;y];
p = polyfit(x,y,2);
y2 = polyval(p,x);
plot(x,y,'b*',x,y2, 'k');
xlabel ('Segment Number')
ylabel ('Tempreture')


%% other functions

function T_res = reservoir_temp (T_top, T_bottom, node_num)
depth = 1632;
dz = 200;
slope = (T_bottom - T_top) / depth;
if node_num*dz < depth
    T_res = slope * node_num*dz + T_top;
else
    T_res = T_bottom;
end
end


function f_tp = friction_tp (alpha_g, alpha_w, rho_mixture, V_mixture)
alpha_o = 1 - alpha_g - alpha_w;
Re = (rho_mixture*V_mixture*0.0762)/(0.001*(alpha_w*1+alpha_o*1));
f_tp = 24/ Re;
end


function Vg = gas_velocity (C0, V_d, mixture_velocity)
Vg = C0 * mixture_velocity + V_d;
end


function Temp = SOR_method (A,B,w, initial_guess, error)
for i = 1: length (A)
    A (i, :) = A (i, :)./ max (abs(A (i, :)));
    B (i, :) = B (i, :)./ max (abs(A (i, :)));
end
X = initial_guess;
cnt = 1;
X_old = X;
while true
    for i = 1:length(A)
        X_P = (B(i) - sum (A(i,:).*X (:)') + A(i,i)* X (i))/A(i,i); 
        X (i) = (1-w)*X (i) + w*X_P; 
    end
    X =X ./ mean(X);
    if cnt == 1 
        error1 = sum (X_old (:));
        error2 = sum(X (:));
        cnt = cnt + 1;
    else
        if abs (error2 - error1) <= error
            Temp = X;
            break
        elseif abs (error2 - error1) > error    
            error1 = sum (X(:));
            cnt = cnt + 1;
        end
    end
end
end

function Uto = depth_heat_coefficient (U_to_top, U_to_end, node_num)
depth = 1632;
slope = (U_to_end - U_to_top)/depth;
dz = 200;
if node_num*dz < depth
    Uto = slope * node_num*depth + U_to_top ;
else
    Uto = U_to_end;
end
end

function U_p = Up (Cp, Tw, T_ref)
U_p = Cp * (Tw - T_ref);
end

function X = rescaling (X)
for i = 1: length (X)
    if mod (i,5) == 0
        X (i) = X(i)*10^(floor(log10(abs(1/X(i))))+1);
    end
end
end

function Q = Q_loss (T_top, T_bottom, U_to_top, U_to_end, T, node_num)
T_res = reservoir_temp (T_top, T_bottom, node_num);
U = depth_heat_coefficient (U_to_top, U_to_end, node_num);
Q = -2*pi*(0.0762/2)*U*(T_res - T);
end
