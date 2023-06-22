function sol = simulation_ODE_model_T(p,tspan)

%options = odeset('RelTol',1e-6,'AbsTol',1e-6,'MaxStep',1e-2,'Events',@events);
options = odeset('RelTol',1e-6,'AbsTol',1e-6,'MaxStep',1e-2); %odeset('NonNegative',1);
sol = ode15s(@ODEs, tspan, p.IC, options);

time = linspace(tspan(1),tspan(end),1000);
%------------------------------------------------------------------------
function dxdt = ODEs(t,x)

%defining the variables
CART_f = x(1); %Free CART concentration
CART_b = x(2); %Bound CART concentration
T = x(3); %Tumor cells concentration
I = x(4); %IFN-gamma concentration
IL = x(5); %Il-6 concentration
Mp = x(6); %Macrophages concentration

%Initialise administration parameters
Infusion=0;%Chemotherapy infusion parameter
Dose=0;%IFN administration parameter

if p.DayAdminChemo <= t && t <= p.DayAdminChemo + (p.Period*p.NumAdminsChemo)%If time is within chemo infusion time
    tchemomod=mod(t-p.DayAdminChemo,p.Period);%Calculate modulo start of infusion time
    %disp('works')
    if p.AdminChemo==0
        Infusion=0;%Infusion=0 since no chemotherapy is administered (p.AdminChemo=0)
    elseif p.AdminChemo==1 && tchemomod < p.DeltaC
        Infusion=p.TotalDose/p.DeltaC;%Administer chemotherapy and calculate dose over infusion time
        %disp('works2')
    end
   tgcsfmod=mod(tchemomod-p.AdminDay,p.IFNPeriod);%Calculate time since chemotherapy admin began
   tGCSFend=p.AdminDay+p.IFNPeriod*p.AdminsIFN;%Calculate time when IFN administrations will end
    if p.AdminDay <= tchemomod && tchemomod <= tGCSFend
        if p.AdminIFN == 0
            Dose=0;%Administer no IFN (p.AdminIFN=0)
        else
            Dose=p.ka*p.F*(p.Dose/p.Vd)*exp(-p.ka*tgcsfmod);%Calculate IFN dose to administer (subcutaneous)
        end
    end
end

%defining the equations
dxdt = [Infusion - p.k_b*(I/p.I0)*(T - CART_b)*(CART_f)^(p.power) + p.k_f*CART_b - p.alpha_T*(I/p.I0)*CART_f*CART_b - (CART_f)/(p.tau_T); %dCART_f/dt
        p.k_b*(I/p.I0)*(T - CART_b)*(CART_f)^(p.power) - p.k_f*CART_b - p.alpha_T*(I/p.I0)*(CART_b)^2 - (CART_b)/(p.tau_T); %dCART_b/dt
        T*p.rho*log((p.T_max)/T) - p.alpha_T*(I/p.I0)*(CART_b)*T;% - p.beta*Mp*T; %dT/dt
        p.r_IFN*CART_b*T + p.P_IFN - p.d_IFN*I; %dIFN/dt
        p.r_IL*(CART_b)*T + (p.r_IL_Mp*Mp)/(Mp/p.Mp0 + p.eta) + p.P_IL - p.d_IL*IL;
        Mp*(p.rho_IL*IL)/(IL + p.epsilon_IL) - p.d_Mp*Mp]; %dMp/dt        
        
end
%------------------------------------------------------------------------
% function [value,isterminal,direction] = events(t,x)
% % Locate the time when tumor passes below a threshold and stop integration to induce tumor eradication.  
% value = x(3) - 0.25;     % detect tumor = 1
% isterminal = 1;   % 1= stop the integration, 0 = continue
% direction = -1;   % negative direction
% end
%------------------------------------------------------------------------
end