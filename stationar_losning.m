 %% Uppgift 1.1 med derivata = 0
 
clear all 
clc
format long
% Fysiskaliska paramtetrar
d=1;          % Längden
rho = 1500;   % Densiteten
lambda = 1.0; % Värmeldningsförmågan
c  = 1000;    % Värmekapaciteten
k = lambda/(c*rho); 
K = 0.50;      % K skall vara av storleken 0.25-0.5 för stabilitet
m=30;         % Antal diskretiseringspunkter
dx=d/(m-1);
dt=K*dx^2/k;


alpha=k*dt/dx^2; % Konstanten framför andra termen i differensekvationen
T_old=ones(m,1)*273; % Temperaturen vid den gamla tidpunkten
T_new=zeros(m,1); % Temperaturen vid den nya tidpunkten
T_old(1) = 263;
T_old(end) = 273;

% Villkor för stationär lösning: Derivatan eller Laplacianen = 0

 for i=2:m-1 % Andra element till näst sista ty första och sista elementet bestäms av randvillkoret
      T_new(i)=T_old(i)+alpha*(T_old(i+1)-2*T_old(i)+T_old(i-1));
 end
 T_new(1) = 263; 
 T_new(end) = 273;
 T_new1 = T_new;
 T_old1 = T_old;
 t = dt;
 
 while all(calcDerivative(T_new1,T_old1,dt) < 1e-14) == false %% Villkor för att hitta stationärlösning
 
    t = t + dt;
    for i=2:m-1
        T_new(i)=T_old(i)+alpha*(T_old(i+1)-2*T_old(i)+T_old(i-1));
    end
    
    % Dirchlets randvillkor
    T_new(1) = 263; 
    T_new(end) = 273;
    
    T_new1 = T_new;
    T_old1 = T_old;
    T_old=T_new;
    
 end
 
 x=linspace(0,d,m);
 plot(x,T_old);
 p = polyfit(x,T_old,1);
 legend
 grid on 
 
 
 %% Uppgift 1.2 med derivata = 0
 
clear all 
clc
format long
% Fysiska paramtetrar
d=1;          % Längden
rho = 1500;   % Densiteten
lambda = 1.0; % Värmeldningsförmågan
c  = 1000;    % Värmekapaciteten
k = lambda/(c*rho); 
s = 100;      % Homogen värmekälla
u = s/(c*rho);% Skalärfält
K = 0.5;      % K skall vara av storleken 0.25-0.5 för stabilitet

m=30;         % Antal diskretiseringspunkter
dx=d/(m-1);
dt=K*dx^2/k;


alpha=k*dt/dx^2; % Konstanten framför andra termen i differensekvationen

T_old=ones(m,1)*273; % Temperaturen vid den gamla tidpunkten
T_new=zeros(m,1); % Temperaturen vid den nya tidpunkten
T_old(1) = 263;
T_old(end) = 273;

% Villkor för stationär lösning: Derivatan eller Laplacianen = 0

 for i=2:m-1 % Andra element till näst sista ty första och sista elementet bestäms av randvillkoret
      T_new(i)=T_old(i)+alpha*(T_old(i+1)-2*T_old(i)+T_old(i-1));
 end
 T_new(1) = 263; 
 T_new(end) = 273;
 T_new1 = T_new;
 T_old1 = T_old;
 t = dt;
 while all(calcDerivative(T_new1,T_old1,dt) < 1e-14) == false %% Villkor för att hitta stationärlösning
 
    t = t + dt;
    for i=2:m-1 % Andra element till näst sista ty första och sista elementet bestäms av randvillkoret
        T_new(i)=T_old(i)+alpha*(T_old(i+1)-2*T_old(i)+T_old(i-1))+dt*u;
    end
    
    % Dirchlets randvillkor
    T_new(1) = 263; 
    T_new(end) = 273;
    
    T_new1 = T_new;
    T_old1 = T_old;
    T_old=T_new;
    
 end
 
 x = linspace(0,d,m);
 p = polyfit(x,T_old,2);
  plot(x,T_old);
 legend
 grid on

 
 
 
 
 
 


 
 
 
 
 
 
 
 
 