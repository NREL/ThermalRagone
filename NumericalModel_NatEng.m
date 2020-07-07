%Numerical model used to quantify the performance of the thermal energy
%stroage device described in the paper titled 'Rate Capability and Ragone
%Plots for Thermal Energy Storage.  Please see read-me file for more
%information.

%Written by: Allison Mahvi
%06/17/2020

tic
%% Operating Conditions
    Control = 'Constant Power (dT)';                                        %Control scheme - Constant Power (dT), Constant Power (m) or Constant Inlet
    T_cutoff = 12 + 273.15;                                                 %Cutoff temperature (for Ragone plots)
    
    if isequal(Control, 'Constant Power (dT)')
        DeltaT_f = 4.5;                                                     %Fluid temperature difference across device (K)
        C_rate = 1;                                                         %Discharge rate (-)
        DisTime = (3600/C_rate);                                            %Discharge time (s)
      	T_PCM_ini = 273.15;                                                 %Initital PCM temperature (K)
        T_f_ini = 273.15;                                                   %Initial Fluid temperture (K)
    elseif isequal(Control, 'Constant Power (m)')
        T_f_inlet = 15 + 273.15;                                            %Fluid inlet temperature (K)
        m_dot_f_max = 0.3;                                                  %Maximum fluid flow rate (kg/s)
        C_rate = 3;                                                         %Discharge rate (-)
        DisTime = 3600/C_rate;                                              %Discharge time (s)
      	T_PCM_ini = 273.15;                                                 %Initital PCM temperature (K)
        T_f_ini = 273.15;                                                   %Initial fluid temperture (K)
    elseif isequal(Control, 'Constant Inlet')
        T_f_inlet = 15 + 273.15;                                            %Specified fluid inlet temperature (K)
        m_dot_f = 2e-2;                                                     %Fluid flow rate (kg/s)
        DisTime = 4000;                                                     %Run simulation for 1 hour (s)
        T_PCM_ini = 273.15;                                                 %Initital PCM temperature (K)
        T_f_ini = 273.15;                                                   %Initial fluid temperture (K)
    else
        disp('Error - Control Scheme Not Coded');                           %Display to user that an error occured     
        return;                                                             %Stop code if control scheme is not properly defined
    end

%% Graphical Outputs
    GeneratePlots = 'yes';                                                   %Generate plots from run?
    
%% Geometry
    L_HX = 0.4572;                                                          %Heat exchanger length (m) - 18 inches total
    W_HX = 0.254;                                                           %Heat exchanger width (m) - 10 inches total
    H_ch = 0.003;                                                           %Channel height (m)
    H_PCM = 0.035;                                                          %Height of PCM (m)
    th_wall = 0.35e-3;                                                       %Wall thickness (m)
    D_h = (4*W_HX*H_ch)/(2*W_HX + 2*H_ch);                                  %Channel hydraulic diameter (m)
    Aspect = H_ch/W_HX;                                                     %Channel aspect ratio (-)
    A_cs = W_HX*H_ch;                                                       %Channel cross-sectional area (m^2)
    rrough = 2.575e-4;                                                      %Relative roughness of channel (-)
    R_contact = 4e-4;                                                       %Contact resistance (K-m^2/W)
    
%% Thermo-Physical Properties
    %Phase Change Composite
    load('h-T_Tetradecane_Graphite.mat');                                   %Load table to find phase change material properties
    k_PCM_ll = 18;                                                          %PCM thermal conductivity in parallel direction (W/mK)
    k_PCM_T = 10;                                                           %PCM thermal conductivity in perpendicular direction (W/mK)
    L = 167977.905255255;                                                   %Latent heat (J/kg)
    rho_PCM = 836;                                                          %PCM density (kg/m^3)
    M_PCM_total = rho_PCM*L_HX*W_HX*H_PCM*2;                                %Total mass of PCM (kg)
    cp_l = 2150;                                                            %Liquid specific heat (J/kg-K)
    cp_s = 2400;                                                            %Solid specific heat (J/kg-K)
    T_bar_pc = 277.75;                                                      %Average phase-change temperature (K)
    dT_glide = 0.5;                                                         %Temperature glide (K)
    T_pc_min = T_bar_pc - (dT_glide/2);                                     %Temperature where phase change initiates (K)
    T_pc_max = T_bar_pc + (dT_glide/2);                                     %Temperature where phase change completes (K)
  	T_ref = 260;                                                            %Reference temperature (K)
    h_ref = 387731.18;                                                      %Reference enthalpy (J/kg)
    
    %Wall
    k_wall = 200;                                                           %Thermal conductivity of tube wall material (W/mK)
    
%% Capacity and Target Heat Transfer Rate
    h_charged = interp1(T_data,h_data,273.15,'linear');                     %Enthalpy where TES is considered fully charged (0C)
    h_discharged = interp1(T_data,h_data,283.15,'linear');                  %Enthalpy where TES is considered fully discharged (10C)
   	h_sat_s = cp_s*(T_pc_min - T_ref) + h_ref;                              %Solid saturation enthalpy of PCM (J/kg)
    h_sat_l = h_sat_s + L;                                                  %Liquid saturation enthalpy of PCM (J/kg)
    dh_dc = h_discharged - h_charged;                                       %Enthalpy change between charge and discharge state (J/kg)
   
    %PCM Capacity
    Cap = (h_discharged - h_charged)*M_PCM_total;                           %Total capacitiy of TES unit cell (J)
 	Q_target = Cap/(3600/C_rate);                                           %Target heat transfer rate (W)
    
%% Discritization
    N_x = 40;                                                               %Number of nodes in fluid direction
    N_y = 20;                                                               %Number of nodes in PCM direction
    dx = L_HX/N_x;                                                          %x-dimension of nodes (m)
    dy = H_PCM/N_y;                                                         %y-dimension of nodes (m)
    x = dx/2:dx:L_HX-dx/2;                                                  %x location of node (m)
    y = dy/2:dy:H_PCM-dy/2;                                                 %y location of node (m)
  	M_PCM = rho_PCM*(dx*W_HX*dy);                                           %Mass of PCM nodes (kg)
  	Conc_f(1:N_x,1) = 10;                                                   %Propylene-Glycol Concentration (%)

%% Time Step
    %Temperature Range
	T_f_crit = 315;                                                         %Highest assumed fluid temperature for critical time step calculation (K)
    T_f_crit_min = T_f_ini;                                               %Lowest assumed fluid temperature for critical time step calculation (K)
    
    %Mass Flow Rate
    if isequal(Control, 'Constant Power (m)')
        m_dot_f = m_dot_f_max;                                              %Worst case scenario for time steps
    end
    if isequal(Control, 'Constant Power (dT)')
        cp_bar_f = PG_Cp([315,Conc_f(1)]')*1000;                            %Specific Heat at High Temperature (J/kg-K)
       	m_dot_f = Q_target/(DeltaT_f*cp_bar_f);                             %Fluid flow rate to calculate dt (High)
    end
    
    %Thermo-Physical Properties to Calculate dt_crit
   	cp_PCM_min = min(cp_l, cp_s);                                           %Minimum specific heat in PCM
  	cp_f_min = PG_Cp([271.15,Conc_f(1)]')*1000;                             %Specific Heat Lowest at Low Temperatures - Assume -2C
    cp_f_max = PG_Cp([T_f_crit,Conc_f(1)]')*1000;                           %Specific Heat Highest at High Temperatures - Assume 315 K
    rho_f_min = PG_rho([T_f_crit,Conc_f(1)]');                              %Fluid Density Lowest at Higher Temperatures - Assume 315 K
 	mu_f_crit = PG_mu([T_f_crit,Conc_f(1)]');                               %Fluid Viscosity
    k_f_crit = PG_k([T_f_crit,Conc_f(1)]');                                 %Fluid thermal conductivity
    Pr_f_crit = PG_Cp([T_f_crit,Conc_f(1)]')*1000*mu_f_crit/k_f_crit;       %Fluid Prandtl number
  	Re_f_crit = 2.*m_dot_f./((W_HX + H_ch).*mu_f_crit);                     %Fluid Reynolds number
  	[Nus_T_crit, Nus_H_crit, f_crit] = ductflow(Re_f_crit,Re_f_crit,...     %Function to Find Single-Phase Heat Transfer Coefficient in a Rectangular Duct
    	Pr_f_crit,L_HX/D_h,Aspect,rrough);
  	htc_f_max = ((Nus_T_crit + Nus_H_crit)/2)*k_f_crit/D_h;                 %Average Fluid Heat Transfer Coefficient in Segment
    
    %Critical Time Steps
   	dt_crit_cent = rho_PCM*cp_PCM_min/2*(k_PCM_ll/dx^2+k_PCM_T/dy^2)^(-1);  %Critical Time Step - PCM Central Node
    dt_crit_edge = rho_PCM*cp_PCM_min*(2*k_PCM_ll/dx^2 + k_PCM_T/dy^2 ...   %Critical Time Step - PCM Node Adjacent to Channel
        +1/(dy*(dy/(2*k_PCM_T)+R_contact+th_wall/k_wall+1/htc_f_max)))^(-1);
    dt_crit_f = rho_f_min*cp_f_min*(m_dot_f*cp_f_max/(dx*H_ch*W_HX) +...    %Critical Time Step - Fluid Nodes
        2/(H_ch*(dy/(2*k_PCM_T)+R_contact+th_wall/k_wall+1/htc_f_max))...
        + 2*k_f_crit/dx^2)^-1;
    dt_crit = min([dt_crit_cent dt_crit_edge 0.8*dt_crit_f]);               %Critical Time Step - Minimum of All Nodes
    
    %Time Step in Simulation
    dt = dt_crit;                                                           %Set Time Step to the Critical Time Step
 
%% Allocate Arrays
  	time = (0:dt:DisTime+dt).';                                             %Time - Allocation
  	timesteps = size(time,1);                                               %Total number of time steps - Allocation
    k_f = nan(N_x,1);                                                     %Fluid thermal conductivity - Allocation
    mu_f = nan(N_x,1);                                                    %Fluid viscosity - Allocation
    Pr_f = nan(N_x,1);                                                    %Fluid Prandtl number - Allocation
    T_bar_f = nan(timesteps,1);                                           %Average fluid temperature at each time step - Allocation
    htc_bar_f = nan(timesteps,1);                                         %Average heat transfer coefficient of fluid - Allocation
    DeltaP_f = nan(timesteps,1);                                          %Average pressure drop across heat exchanger - Allocation
    q_f_PCM = nan(N_x,timesteps);                                         %Boundary heat transfer rate - Allocation
    dhdt = nan(N_x,N_y,timesteps);                                        %Enthalpy change of PCM with time - Allocation
    h = nan(N_x,N_y,timesteps);                                           %Enthalpy of PCM - Allocation
    T = nan(N_x,N_y,timesteps);                                           %Temperature of PCM - Allocation
    T_f = nan(N_x,timesteps);                                             %Temperature of fluid - Allocation
    T_f_in = nan(timesteps,1);                                            %Inlet fluid temperature - Allocation
    q_f_PCM_total = nan(timesteps,1);                                     %Total heat transfer from fluid to PCM - Allocation
    q_f_total = nan(timesteps,1);                                         %Total fLuid heat transfer (inlet/outlet) - Allocation
    h_bar = nan(timesteps,1);                                             %Average PCM enthalpy - Allocation
    SOC = nan(timesteps,1);                                               %State of charge - Allocation
    q_flux_surf = nan(N_x,timesteps);                                     %Surface heat flux - Allocation
    dTdt_f = nan(N_x, timesteps-1);                                       %Fluid temperature change with time - Allocation
    q_f_stored = nan(timesteps-1, 1);                                     %Fluid heat storage - Allocation
	m_f_track = nan(timesteps,1);                                         %Track flow rate with time - Allocation
    Phase = nan(N_x,N_y);                                                 %Phase of each node - Allocation
    
    %Resistance/Area Calculations
   	Nodes_active = nan(N_x,1);                                              %Active x-locations - Allocation
    Nodes_liquid = nan(N_x,1);                                              %x-locations with a liquid layer - Allocation
    AreaRatio = nan(timesteps,1);                                           %Active area - Allocation
    Area_liquid = nan(timesteps,1);                                         %Area with liquid layer - Allocation
    L_front = nan(N_x,1);                                                   %Distance between channel and liquid front - Allocation
    R_fluid = nan(timesteps,1);                                             %Fluid resistance per unit area - Allocation
    R_liquid = nan(timesteps,1);                                            %Liquid PCM resistance - Allocation
    R_wall = nan(timesteps,1);                                              %Wall resistance - Allocation
    R_cont = nan(timesteps,1);                                              %Contact resistance - Allocation
    
%% Set Values for First Time Step
    %Thermophysical Properties
    T(1:N_x,1:N_y,1) = T_PCM_ini;                                           %Initial condition for PCM nodes (K)
    h(:,:,1) = interp1(T_data,h_data,T(:,:,1),'linear');                    %Inital PCM enthalpy (J/kg)     
 	T_f(1:N_x,1) = T_f_ini;                                                 %Initial condition for fluid nodes (K)
	T_bar_f(1) = mean(T_f(:,1));                                            %Average fluid temperature (K)

    %PCM Phase
    Phase(h(:,:,1) < h_sat_s) = 0;                                          %Solid nodes
    Phase(h(:,:,1) > h_sat_l) = 1;                                          %Liquid nodes
    Phase(h(:,:,1) > h_sat_s & h(:,:,1) < h_sat_l) = 0.5;                   %Mushy nodes
	Phase_wall = Phase(:,1);                                                %Phase at wall
	Phase_boundary = Phase(:,N_y);                                          %Phase at far boundary
    
    %Inlet Conditions
    if isequal(Control, 'Constant Power (dT)')
        T_f_in(1) = T_f(N_x,1) + DeltaT_f;                                  %Initial inlet fluid temperature (K)
    end
    if isequal(Control,'Constant Power (m)') || isequal(Control,'Constant Inlet')
        T_f_in(1:timesteps) = T_f_inlet;                                    %Initial inlet fluid temperature (K)
    end
    
%% Calculate mass flow rates
    if isequal(Control, 'Constant Power (dT)')
        cp_bar_f = PG_Cp([T_bar_f(1),Conc_f(1)]')*1000;                     %Average specific heat of the fluid
        m_dot_f = Q_target/(DeltaT_f*cp_bar_f);                             %Correct initial fluid temperature (K)
        m_f_track(1) = m_dot_f;                                             %Track flow rate with time (kg/s)

    elseif isequal(Control, 'Constant Power (m)')
        cp_bar_f = PG_Cp([T_bar_hot(1),Conc_f(1)]')*1000;                   %Average specific heat of the fluid
        m_dot_f = Q_target/(cp_bar_f*(T_f_in(1) - T_f(N_x,1)));             %Target fluid mass flow rate (kg/s)
        m_f_track(1) = m_dot_f;                                             %Track flow rate with time (kg/s)
    end
    
%% Output Information for Users
    F_disp = 0.1;                                                           %Fraction of total time when elapsed time is displayed                                                                                      
    dispinc = 0.1;                                                          %Increments that elapsed time is displayed     
    
%% Loop Through Time
for t = 1:1:timesteps-1
     %% Provide Information to User About Progress
        if t == ceil(timesteps*F_disp)
            disp([num2str(F_disp*100),'% Complete - Time Elapsed '...       %Information for user
                ,num2str(toc),' sec']);
            F_disp = F_disp + dispinc;                                      %Next fraction that will be displayed
        end
  
    %% State of Charge
        h_bar(t) = mean(h(:,:,t),'all');                                    %Average PCM enthalpy (J/kg)
        SOC(t) = 100-((h_bar(t)-h_charged)/(h_discharged-h_charged))*100;   %State of charge (%)
        
    %% Fluid Thermophysical Properties
        %Average Fluid Properties (Change less than 1%)
        T_bar_f(t) = mean(T_f(:,t));                                        %Average fluid temperature (K)
        cp_bar_f = PG_Cp([T_bar_f(t),Conc_f(1)]')*1000;                     %Average specific heat of the fluid
       	rho_bar_f = PG_rho([T_bar_f(t),Conc_f(1)]');                        %Fluid density (averaged over x)
        M_bar_f = rho_bar_f*dx*H_ch*W_HX;                                   %Average mass of fluid node (averaged over x)
        k_bar_f = PG_k([T_bar_f(t),Conc_f(1)]');                            %Average fluid thermal conductivity
  
      	%Local Hot Fluid Properties
        k_f = PG_k([T_f(:,t), Conc_f]')';                                   %Fluid thermal conductivity
        mu_f = PG_mu([T_f(:,t), Conc_f]')';                                 %Fluid viscosity
        Pr_f = cp_bar_f.*mu_f./k_f;                                         %Fluid Prandtl number
        Re_f = 2.*m_dot_f./((W_HX + H_ch).*mu_f);                           %Fluid Reynolds number
        
        %Averaged Flow Properties
        Re_f_bar = mean(Re_f);                                              %Average fluid reynolds number
        u_m = m_dot_f/(rho_bar_f*W_HX*H_ch);                                %Average velocity

    %% Heat Transfer Coefficient and Pressure Drop
        %Local
        if not(isequal(Control, 'No Flow'))
            [Nus_T, Nus_H, f] = ductflow(Re_f_bar,Re_f(:),Pr_f(:),...       %Function to find single-phase heat transfer coefficient in a rectangular duct
                L_HX/D_h,Aspect,rrough);
            Nus_bar = (Nus_T + Nus_H)./2;                                   %Average fluid Nusselt number (-)
            htc_f = Nus_bar.*k_f./D_h;                                      %Average fluid heat transfer coefficient in segment (W/m^2K)
            dPdx = f.*(1/D_h).*(1/2).*rho_bar_f.*u_m^2;                     %Pressure gradient (Pa/m)
        else
            htc_f(1:N_x,1) = 5;                                             %Natural convective heat transfer coefficient (Guess - W/m^2K)
        end
        
        %Average
        htc_bar_f(t) = mean(htc_f(1:N_x));                                  %Average Fluid Heat Transfer Coefficient across HX
        DeltaP_f(t) = mean(dPdx(1:N_x))*L_HX;                               %Average pressure drop across heat exchanger
    
    %% Fluid/PCM Interactions
        q_f_PCM(1:N_x,t) = dx.*((T_f(1:N_x,t) - T(1:N_x,1,t))./...          %Heat transfer rate in each segement from fluid to PCM
            (dy/(2*k_PCM_T*W_HX) + (R_contact/W_HX)+...
            (th_wall/(k_wall*W_HX)) + 1./(htc_f(1:N_x).*W_HX)));
        q_f_PCM_total(t) = sum(q_f_PCM(1:N_x,t))*2;                         %Total heat transfered from fluid to PCM in each time step
        q_f_total(t) = m_dot_f*cp_bar_f*(T_f_in(t) - T_f(N_x,t));           %Total fluid heat transfer rate (using inlet and outlet temperatures)
        q_flux_surf(:,t) = q_f_PCM(:,t)./(dx*W_HX);                         %Surface heat flux
        
    %% Energy Balances
        %PCM - Insulated Corner Nodes (Top Left and Right)
        dhdt(1,N_y,t)=1/M_PCM .*(k_PCM_ll*dy*W_HX/dx.*(T(2,N_y,t)-T(1,N_y,t)) + k_PCM_T*dx*W_HX/dy.*(T(1,N_y-1,t)-T(1,N_y,t)));
        dhdt(N_x,N_y,t) = 1/M_PCM *(k_PCM_ll*dy*W_HX/dx*(T(N_x-1,N_y,t)-T(N_x,N_y,t)) + k_PCM_T*dx*W_HX/dy*(T(N_x,N_y-1,t)-T(N_x,N_y,t)));
    
        %PCM - Corner Nodes Adjacent to Channels (Bottom Left and Right)
        dhdt(1,1,t) = 1/M_PCM .*(k_PCM_ll*dy*W_HX/dx.*(T(2,1,t)-T(1,1,t)) + k_PCM_T*dx*W_HX/dy.*(T(1,2,t)-T(1,1,t)) + q_f_PCM(1,t));
        dhdt(N_x,1,t) = 1/M_PCM *(k_PCM_ll*dy*W_HX/dx*(T(N_x-1,1,t)-T(N_x,1,t)) + k_PCM_T*dx*W_HX/dy*(T(N_x,2,t)-T(N_x,1,t)) + q_f_PCM(N_x,t));
    
        %PCM - Top Insulated Edge
        dhdt(2:N_x-1,N_y,t) = 1/M_PCM *(k_PCM_ll*dy*W_HX/dx*(T(1:N_x-2,N_y,t)-T(2:N_x-1,N_y,t)) + k_PCM_ll*dy*W_HX/dx*(T(3:N_x,N_y,t)-T(2:N_x-1,N_y,t)) + k_PCM_T*dx*W_HX/dy*(T(2:N_x-1,N_y-1,t)-T(2:N_x-1,N_y,t)));

        %PCM - Bottom Edge Adjacent to Fluid Channel
        dhdt(2:N_x-1,1,t) = 1/M_PCM *(k_PCM_ll*dy*W_HX/dx*(T(1:N_x-2,1,t)-T(2:N_x-1,1,t)) + k_PCM_ll*dy*W_HX/dx*(T(3:N_x,1,t)-T(2:N_x-1,1,t)) + k_PCM_T*dx*W_HX/dy*(T(2:N_x-1,2,t)-T(2:N_x-1,1,t)) + q_f_PCM(2:N_x-1,t));
    
        %PCM - Left and Right Edges
        dhdt(1,2:N_y-1,t) = 1/M_PCM *(k_PCM_ll*dy*W_HX/dx*(T(2,2:N_y-1,t)-T(1,2:N_y-1,t)) + k_PCM_T*dx*W_HX/dy*(T(1,3:N_y,t)-T(1,2:N_y-1,t)) + k_PCM_T*dx*W_HX/dy*(T(1,1:N_y-2,t)-T(1,2:N_y-1,t)));
        dhdt(N_x,2:N_y-1,t) = 1/M_PCM *(k_PCM_ll*dy*W_HX/dx*(T(N_x-1,2:N_y-1,t)-T(N_x,2:N_y-1,t)) + k_PCM_T*dx*W_HX/dy*(T(N_x,3:N_y,t)-T(N_x,2:N_y-1,t)) + k_PCM_T*dx*W_HX/dy*(T(N_x,1:N_y-2,t)-T(N_x,2:N_y-1,t)));

        %PCM - Central Nodes
        dhdt(2:N_x-1,2:N_y-1,t) = 1/M_PCM *(k_PCM_ll*dy*W_HX/dx*(T(1:N_x-2,2:N_y-1,t) - T(2:N_x-1,2:N_y-1,t)) + k_PCM_ll*dy*W_HX/dx*(T(3:N_x,2:N_y-1,t) - T(2:N_x-1,2:N_y-1,t)) + k_PCM_T*dx*W_HX/dy*(T(2:N_x-1,3:N_y,t) - T(2:N_x-1,2:N_y-1,t)) + k_PCM_T*dx*W_HX/dy*(T(2:N_x-1,1:N_y-2,t) - T(2:N_x-1,2:N_y-1,t)));
        
        %Fluid Energy Balance (With Fluid Conduction)
        dTdt_f(1,t) = 1/(cp_bar_f*M_bar_f)*(m_dot_f*cp_bar_f*(T_f_in(t)-T_f(1,t)) - 2*q_f_PCM(1,t) + k_bar_f*H_ch*W_HX/dx*(T_f(2,t) - T_f(1,t)));
        dTdt_f(2:N_x-1,t) = 1/(cp_bar_f*M_bar_f)*(m_dot_f*cp_bar_f*(T_f(1:N_x-2,t)-T_f(2:N_x-1,t)) - 2*q_f_PCM(2:N_x-1,t) + k_bar_f*H_ch*W_HX/dx*(T_f(1:N_x-2,t)-T_f(2:N_x-1,t)) + k_bar_f*H_ch*W_HX/dx*(T_f(3:N_x,t) - T_f(2:N_x-1,t)));
        dTdt_f(N_x,t) = 1/(cp_bar_f*M_bar_f)*(m_dot_f*cp_bar_f*(T_f(N_x-1,t)-T_f(N_x,t)) - 2*q_f_PCM(N_x,t) + k_bar_f*H_ch*W_HX/dx*(T_f(N_x-1,t) - T_f(N_x,t)));

        q_f_stored(t) = (cp_bar_f*M_bar_f)*sum(dTdt_f(1:N_x,t));     
 
    %% Liquid PCM Resistance and Active Area
        %Liquid and Active Cross Sectional Areas
        Nodes_liquid(Phase_wall < 1) = 0;                                   %x-locations with a liquid layer
        Nodes_liquid(Phase_wall == 1) = 1;                                  %x-locations with a liquid layer
        Nodes_active(Phase_boundary < 1) = 1;                               %x-nodes that are active
        Nodes_active(Phase_boundary == 1) = 0;                              %x-nodes that are inactive
        N_x_active = sum(Nodes_active(:));                                  %Total number of active nodes
        N_x_liquid = sum(Nodes_liquid(:));                                  %Total number of inactive nodes
        Area = N_x_active*dx*W_HX;                                          %Active area (m^2)
        AreaRatio(t) = (Area/(L_HX*W_HX))*100;                              %Active area ratio (%)
        Area_liquid(t) = N_x_liquid*dx*W_HX;                                %Area with a liquid layer (m^2)
        
        %Location of Phase Front
      	N_liquid = zeros(N_y,1);                                            %Boolean - liquid node?
        N_y_front = zeros(N_x,1);                                           %y-node of liquid front
        for i = 1:N_x
            N_liquid(Phase(i,:) < 1) = 0;                                   %Solid or mushy nodes
            N_liquid(Phase(i,:) == 1) = 1;                                  %Liquid nodes
            N_y_front(i) = sum(N_liquid);                                   %y-node of liquid front
        end

        %Calculate Resistances
       	R_fluid(t) = 1/htc_bar_f(t);                                        %Fluid resistance per unit area
      	R_PCM_max = H_PCM/k_PCM_T;                                          %Maximum PCM resistance
        R_wall(t) = th_wall/k_wall;                                         %Wall resistance (constant)
        R_cont(t) = R_contact;                                              %Contact resistance (constant)
        if Area_liquid(t) == 0
            R_liquid(t) = 0;                                                %Liquid resistance per unit area (if no liquid)
        else
       	if Area_liquid(t)/(W_HX*L_HX) < 0.9999  
          	L_front(1:N_x_liquid) = N_y_front(1:N_x_liquid)*dy;             %Length of liquid region (m)
            L_bar_front = mean(L_front(1:N_x_liquid));                      %Average length of liquid region (m)
            R_liquid(t) = L_bar_front/k_PCM_T;                              %Liquid resistance
        else
       	    L_front(1:N_x) = N_y_front(1:N_x)*dy;                           %Length of liquid region (m)
        	L_bar_front = mean(L_front(N_x-N_x_active+1:N_x));              %Average length of liquid region (m)
            R_liquid(t) = L_bar_front/k_PCM_T;                              %Liquid resistance
       	end
        end
    
    %% Integration over Time
        h(1:N_x,1:N_y,t+1) = h(1:N_x,1:N_y,t) + dhdt(1:N_x,1:N_y,t).*dt;    %PCM enthalpy at next time step
        T(:,:,t+1) = interp1(h_data,T_data,h(:,:,t+1),'linear');            %PCM temperature at next time step
        T_f(1:N_x,t+1) = T_f(1:N_x,t) + dTdt_f(1:N_x,t).*dt;                %Fluid temperature at next time step
        
    	Phase(h(:,:,t+1) < h_sat_s) = 0;                                    %Solid nodes
        Phase(h(:,:,t+1) > h_sat_l) = 1;                                    %Liquid nodes
        Phase(h(:,:,t+1) > h_sat_s & h(:,:,t+1) < h_sat_l) = 0.5;           %Mushy nodes
        Phase_wall = Phase(:,1);                                            %Phase at wall
        Phase_boundary = Phase(:,N_y);                                      %Phase at far boundary

        if isequal(Control, 'Constant Power (dT)')
            T_f_in(t+1) = T_f(N_x,t+1) + DeltaT_f;                          %Inlet fluid temperature at next time step (constant power condition)
            T_bar_f(t+1) = mean(T_f(:,t+1));                                %Average fluid temperature
            cp_bar_f = PG_Cp([T_bar_f(t+1),Conc_f(1)]')*1000;               %Average specific heat of the fluid
            m_dot_f = Q_target/(DeltaT_f*cp_bar_f);                         %Update Mass Flow Rate to Compensate for Changes in cp
            m_f_track(t+1) = m_dot_f;                                       %Track Flow Rate with Time in This Control Case (kg/s)
            
        elseif isequal(Control, 'Constant Power (m)')
        	T_bar_f(t+1) = mean(T_f(:,t+1));                               	%Average fluid temperature
            cp_bar_f = PG_Cp([T_bar_f(t+1),Conc_f(1)]')*1000;               %Average specific heat of the fluid
        	m_dot_f = Q_target/(cp_bar_f*(T_f_in(t) - T_f(N_x,t)));         %Target fluid mass flow rate (kg/s)
            if m_dot_f>m_dot_f_max
                m_dot_f = m_dot_f_max;                                      %Sets flow rate to limit if maximum is reached (kg/s)
            end
            m_f_track(t+1) = m_dot_f;                                       %Track flow rate with time (kg/s)

        elseif isequal(Control, 'Constant Inlet')
            T_bar_f(t+1) = mean(T_f(:,t+1));                                %Average fluid temperature
        end

end
    
 %% Final Global Parameters
    t = timesteps;                                                          %Final Time Evaluated (s)
  	h_bar(t) = mean(h(:,:,t),'all');                                        %Average PCM Enthalpy (J/kg)
    DeltaP_bar = mean(DeltaP_f,'all');                                      %Average fluid pressure drop (Pa)
   	SOC(t) = 100-((h_bar(t)-h_charged)/(h_discharged-h_charged))*100;       %State of charge (%)
  	q_f_PCM(1:N_x,t) = dx.*((T_f(1:N_x,t) - T(1:N_x,1,t))./...              %Heat transfer rate in each segement from fluid to PCM
        (dy/(2*k_PCM_T*W_HX) + (R_contact/W_HX)+...
        (th_wall/(k_wall*W_HX)) + 1./(htc_f(1:N_x).*W_HX)));
   	q_f_PCM_total(t) = sum(q_f_PCM(1:N_x,t))*2;                             %Total heat transfered from fluid to PCM in each time step (W)
	q_f_total(t) = m_dot_f*cp_bar_f*(T_f_in(t) - T_f(N_x,t));               %Total fluid heat transfer rate (W)
    RunTime = toc;                                                          %Time to run simulation (s)
    disp(['Simulation Complete - Total Run Time = ',num2str(toc),' sec']);  %User information
 
 %% Specific Power and Energy
    %Cutoff Time
    if T_f(N_x, timesteps) < T_cutoff
        t_cutoff = time(timesteps);                                         %Cutoff time (s)
        Index_cutoff = timesteps;                                           %Index of cutoff (-)
        SOC_cutoff = SOC(timesteps);                                        %Cutoff state of charge (%)
 	else
        intersect = abs(T_f(N_x,:)- T_cutoff);                              %Makes minimum at cutoff time
        [Int_max,Index_cutoff] = min(intersect);                            %Finds cutoff index (-)
        t_cutoff = time(Index_cutoff);                                      %Cutoff time (s)
        SOC_cutoff = SOC(Index_cutoff);                                     %Cutoff state of charge (%)
    end
    
    %Specific Power
        Pow = mean(q_f_PCM_total(1:Index_cutoff));                          %Average power (W-hr)
        SP = mean(q_f_PCM_total(1:Index_cutoff))/M_PCM_total;               %Gravimetric specific power (PCM - W-hr/kg)
        SP_f = mean(q_f_total(1:Index_cutoff))/M_PCM_total;                 %Gravimetric specific power (fluid - W-hr/kg)
        SP_vol = SP*rho_PCM;                                                %Volumetric specific power (PCM - W-hr/m^3)
        SP_vol_f = SP_f*rho_PCM;                                            %Volumetric specific power (fluid - W-hr/m^3)
        
    %Specific Energy
        E = q_f_PCM_total(1:Index_cutoff).*dt;                              %Energy extracted (W)
        E_f = q_f_total(1:Index_cutoff).*dt;                                %Fluid energy dissipated (W)
        SE = (sum(E)/3600)/M_PCM_total;                                     %Gravimetric specific energy (PCM - W/kg)
        SE_f = (sum(E_f)/3600)/M_PCM_total;                                 %Gravimetric specific energy (fluid - W/kg)
        SE_vol = SE*rho_PCM;                                                %Volumetric specific energy (PCM - W/m^3)
        
%% Generate Plots
    pInc = ceil(3/dt);                                                      %Plot every 3 data points
    T_C = T - 273.15;                                                       %PCM temperatures in C
    T_f_C = T_f - 273.15;                                                   %Fluid temperatures in C
    T_cutoff_C = T_cutoff - 273.15;                                         %Cutoff temepratures in C

    if isequal(GeneratePlots,'yes')
        %General Plot Information
            colors = [60/255 133/255 244/255; 225/255 68/255 55/255;        %Default color scheme
                244/255 160/255 0; 15/255 157/255 88/255; 
                144/255 103/255 167/255; 128/255 133/255 133/255];
           	set(0, 'DefaultLineLineWidth', 2);                              %Set default line width
            set(groot, 'defaultAxesColorOrder', colors);                    %Set default color scheme
            
     	%Rate Capability Plot
        	figure (1);
            plot(SOC(1:pInc:timesteps), T_f_C(N_x,1:pInc:timesteps))
            hold on;
            %plot([0 100], [T_cutoff_C T_cutoff_C],'k--');
            %text(98, T_cutoff_C+0.5, 'Cutoff Temperature');
            ax = gca;
         	ax.FontSize = 12;
            xlim([0 100]);
            set ( gca, 'xdir', 'reverse' )
            xlabel('State of charge (%)', 'FontSize', 14, 'FontWeight', 'bold');
            ylim([0 18]);
            ylabel('Fluid outlet temperaure (°C)', 'FontSize', 14, 'FontWeight', 'bold');
            
      	%Ragone Plot
            figure (2);
            PCMmarker = loglog(SE, SP,'o');
            set(PCMmarker, 'markerfacecolor', get(PCMmarker, 'color'));
           	ax = gca;
         	ax.FontSize = 12;
            hold on
            xlim([10 60]);
            xticks([10 20 30 40 50 60 70]);
           	xlabel('Specific energy (W-hr/kg)', 'FontSize', 14, 'FontWeight', 'bold');
            ylim([8 160]);
            yticks([5 10 20 40 80 120 160 200]);
            ylabel('Specific power (W/kg)', 'FontSize', 14, 'FontWeight', 'bold');
    end
        