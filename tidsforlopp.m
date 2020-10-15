 %% Ritar ut alla tidsförlopp
 
clear all 
clc
clf
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

cell_data = cell(3888,2); % Första kolumn är temperatur och andra tid
j = 0;
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
    j = j + 1;
    t = t + dt;
 
    cell_data{j,1} = T_old1;
    cell_data{j,2} = t;
    
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
 
 grid on 
 hold on 
 x=linspace(0,d,m);
 plot(x,cell_data{1,1}); 

 for i = 1:16
    time(i,1) = cell_data{i*25,2};
    plot(x,cell_data{i*25,1}); 
 end

 
 days = strings(16,1);
 hours = strings(16,1);
 minutes = strings(16,1);
 seconds = strings(16,1);
 
 time(:,1) = sort(round(time(:,1),0));
 days(:,1) = num2str(floor(time(:,1)/(24*3600)));
 hours(:,1) = num2str(round(mod(time(:,1),24*3600)/3600,0));
 minutes(:,1) = num2str(round(mod(mod(time(:,1),24*3600),3600)/60,0));
 seconds(:,1) = num2str(round(mod(mod(mod(time(:,1),24*3600),3600),60),0));



concat = 't = ' + days + ' dagar ' + hours + ' h ' + minutes + ' min ' + seconds + ' s';
legend(concat(:))

 
 %% Ritar ut alla tidsförlopp
 
clear all 
clf
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
 j = 0;
 cell_data = cell(3890,1);
 
 
 while all(calcDerivative(T_new1,T_old1,dt) < 1e-14) == false %% Villkor för att hitta stationärlösning
    j = j+1;
    t = t + dt;
    
    cell_data{j,1} = T_old1;
    cell_data{j,2} = t;
    
    
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
 
 grid on
 hold on
 x=linspace(0,d,m);

 for i = 1:12
     time(i,1) = cell_data{i*25,2};
    plot(x,cell_data{i*25,1}); 
 end
 k = 12;
 days = strings(k,1);
 hours = strings(k,1);
 minutes = strings(k,1);
 seconds = strings(k,1);
 
 time(:,1) = sort(round(time(:,1),0));
 days(:,1) = num2str(floor(time(:,1)/(24*3600)));
 hours(:,1) = num2str(round(mod(time(:,1),24*3600)/3600,0));
 minutes(:,1) = num2str(round(mod(mod(time(:,1),24*3600),3600)/60,0));
 seconds(:,1) = num2str(round(mod(mod(mod(time(:,1),24*3600),3600),60),0));



concat = 't = ' + days + ' dagar ' + hours + ' h ' + minutes + ' min ' + seconds + ' s';
legend(concat(:))

 

 
 
 
 
 
 
 
 
 
 
 
 
 
 