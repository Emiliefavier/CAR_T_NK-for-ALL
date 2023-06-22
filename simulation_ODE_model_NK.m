function sol = simulation_ODE_model_NK(p,tspan)

options = odeset('RelTol',1e-6,'AbsTol',1e-6,'MaxStep',1e-2);
sol = ode15s(@ODEs, tspan, p.IC, options);

time = linspace(tspan(1),tspan(end),1000);
%------------------------------------------------------------------------
function dxdt = ODEs(t,x)

%defining the variables
CARNK_f = x(1); %Free CARNK concentration
CARNK_b = x(2); %Bound CARNK concentration
T = x(3); %Tumor cells concentration
I = x(4); %IFN-gamma concentration
IL = x(5); %Il-6 concentration
G = x(6); %GM-CSF concentration
Mc = x(7); %Monocytes concentration
Mp = x(8); %Macrophages concentration

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
        if p.AdminGCSF == 0
            Dose=0;%Administer no IFN (p.AdminIFN=0)
        else
            Dose=p.ka*p.F*(p.Dose/p.Vd)*exp(-p.ka*tgcsfmod);%Calculate IFN dose to administer (subcutaneous)
        end
    end
end

%defining the equations
dxdt = [Infusion- p.k_b*(I/p.I0)*(T - CARNK_b)*(CARNK_f)^(p.power) + p.k_f*CARNK_b - p.alpha_NK*(I/p.I0)*CARNK_f*CARNK_b - (CARNK_f)/(p.tau_NK); %dCARNK_f/dt
         p.k_b*(I/p.I0)*(T - CARNK_b)*(CARNK_f)^(p.power) - p.k_f*CARNK_b - p.alpha_NK*(I/p.I0)*(CARNK_b)^2 - (CARNK_b)/(p.tau_NK); %dCARNK_b/dt
         T*p.rho*log((p.T_max)/T) - p.alpha_NK*(I/p.I0)*(CARNK_b)*T - p.beta*Mp*T; %dT/dt
         p.r_IFN*CARNK_b*T + p.P_IFN - p.d_IFN*I; %dIFN/dt
         (p.r_IL_Mp*Mp)/(Mp + p.eta)+ p.P_IL - p.d_IL*IL; %dIL/dt
         p.r_G*CARNK_b*T + p.P_G - p.d_G*G; %dG/dt
         (p.Mc_prod*(p.Mc_max - p.Mc_prod)*((G^(p.h_Mc))/(G^(p.h_Mc) + p.epsilon_Mc^(p.h_Mc))))*p.Mc_I - (p.rho_Mcp*G^(p.h_Mc)*Mc)/(G^(p.h_Mp) + p.epsilon_Mp^(p.h_Mp)) + (p.p_M*(T/5000)*Mc)/((T/5000) + p.e_M) - p.d_Mc*Mc; %dMc/dt
         (p.rho_Mcp*G^(p.h_Mc)*Mc)/(G^(p.h_Mp) + p.epsilon_Mp^(p.h_Mp)) - p.d_Mp*Mp];%dMp/dt

end
end
