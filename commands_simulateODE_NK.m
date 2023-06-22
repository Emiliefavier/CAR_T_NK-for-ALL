%CAR NK
%Set initial conditions
p.CARNKf0 = 0;
p.CARNKb0 = 0;
p.T0 = 1; %2 for event;
p.I0 = 10.2; %baseline concentration
p.IL0 = 1.1;
p.Mc0 = 0; 
p.Mp0 = 0;
p.IC = [p.CARNKf0,p.CARNKb0,p.T0,p.I0,p.IL0,p.G0,p.Mc0,p.Mp0];
%Set parameters
p.k_b = 1.3*10^(-3)*86400*10^6; 
p.k_f = 24*86400;
p.power = 1;
p.tau_NK = 10;
p.alpha_NK = 10^(0);
p.rho = 1/45;
p.T_max = 10^3;
p.beta = 10^(-11)*10^9;
p.r_IL_Mp = 1872;
p.d_IL = 16.6;
p.eta = 3.6*10^(-5);
p.P_IL = p.d_IL*p.IL0 - (p.r_IL_Mp*p.Mp0)/(p.Mp0 + p.eta);
p.r_IFN = 0.0315*10^(-6);
p.d_IFN = 1.51;
p.P_IFN = p.d_IFN*p.I0 - p.r_IFN*p.CARNKb0*p.T0;
p.r_G = 5; %assumed
p.d_G = 11.74;
p.P_G = p.d_G*p.G0 - p.r_G*p.CARNKb0*p.T0;
p.Mc_max = 11.55;
p.Mc_I = 0.0023;
p.h_Mc = 1.67;
p.h_Mp = 2.03;
p.epsilon_Mc = 57.2;
p.epsilon_Mp = 0.027;
p.rho_Mcp = 1.7;
p.d_Mp = 0.3;
p.d_Mc = 0.76;
p.Mc_prod= 0.5;
p.d_Mp =0.3;
%------------------------------------------------------------------------
%The next two if statements administer the doses of chemotherapy and cytokine over a given infusion time

%Chemo related parameters 

%Administer CAR (1=yes, 0=no)
p.AdminChemo=1; 
p.AdminGCSF=0;

%Set doses
p.BM=70; %(kg)
p.DoseChemo=10^(-1);%(10^9 cells/kg)
p.TotalDose= p.DoseChemo*p.BM; %(10^9 cells)

%Set regimen
p.Period=16; %time period between the chemos
p.DayAdminChemo=0; %day of chemo administration
p.DeltaC=1/15; %infusion time
p.NumAdminsChemo=1; %number of chemo administrations

%Cytokine related parameters
p.AdminIFN=1; %Administer cytokine(0,=no, 1=yes)
p.AdminDay=0; %day of cytokine administration
p.Dose=3440; %dose
p.IFNPeriod=1; %time period between 2 cytokines administrations
p.AdminsIFN=1; %number of cytokine admnistration

%Dose-Dep SC parameters
p.F=0.89; %bioavailibility of IFN (89%)
p.ka=1/2.5; %absorption constant
p.Vd=5000; %distribution volume
p.ivt0=0; %Time IV infusion begins (in days)
p.tinf=1; %Duration of IV infusion (in days))
%------------------------------------------------------------------------
%Set parameters for the ODE solver
tspan = [0 30]; %time span to simulate over

%Call ODE solver function
sol = simulation_ODE_model_NK(p,tspan);

time = tspan(1):0.01:tspan(2);
sol_mesh = deval(sol,time);

CARNK_f = sol_mesh(1,:);
CARNK_b = sol_mesh(2,:);
T = sol_mesh(3,:);
I = sol_mesh(4,:);
IL = sol_mesh(5,:);
G = sol_mesh(6,:);
Mc = sol_mesh(7,:).*5000;
Mp = sol_mesh(8,:).*5000;
%------------------------------------------------------------------------
%%Plot result 
figure(1)
hold on %to plot multiple plots on the same figure
plot(time,CARNK_f,'Color','cyan','LineWidth',2)
plot(time,CARNK_b,'Color','blue','LineWidth',2)
plot(time,T,'Color','red','LineWidth',2)
xlabel('Time (days)')
ylabel('Number of cells (x10^9)')
set(gca,'FontSize',16)
legend('CAR NK free','CAR NK bound','Tumour cells',Location='best')

figure(2)
hold on
plot(time,I,'Color','yellow','LineWidth',2)
plot(time,IL,'Color','green','LineWidth',2)
plot(time,G,'Color','#7E2F8E','lineWidth',2)
xlabel('Time (days)')
ylabel('Concentration (pg/mL)')
legend('IL-6 concentration',Location='best')
legend('IFN-gamma concentration','IL-6 concentration','GM-CSF concentration',Location='best')

figure(3)
hold on
plot(time,Mp,'Color','magenta')
plot(time,Mc,'Color','#A2142F','LineWidth',2)
xlabel('Time (days)')
ylabel('Number of cells (x10^9 cells)')
legend('Macrophages',Location='best')
