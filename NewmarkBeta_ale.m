%% Incializar
clear variables
close all
clc

%% Datos
m = 0.1; %lb-sec2/in
zeta = 0;
k = 40; %lb/in
syms t
p(t) = 10*cos(10*t); % Excitación (contínua pero se discretiza después)
wn = 20; %rad/sec
Ttot = 0.2; %sec
c = 0;  % Amortiguamiento 
dt = [0.02, 0.005]; % paso

gamma = 1/2;
disp('Están ordenadas para:');
disp('beta = 1/4 dt = 0.02');
disp('beta = 1/4 dt = 0.005');
disp('beta = 1/6 dt = 0.02');
disp('beta = 1/6 dt = 0.005');
disp('No supe ordenarlas con las funciones dentro del espacio/ambiente de figure');
for z = 1:2 % Para recorrer los beta
    if z == 1
        beta = 1/4; %Aceleración promedio
    elseif z == 2
        beta = 1/6; %Aceleración lineal
    end
    for j = 1:length(dt) % Recorriendo los tiempos pasos
        dtj = dt(1,j);
        N = Ttot/dtj; % particiones
        u = zeros(N+1,1);
        v = zeros(N+1,1);
        ac = zeros(N+1,1);
        pu = zeros(N+1,1);
        pv = zeros(N+1,1);
        times = zeros(N+1,1);
        ac(1,1) = (p(0)-c*v(1,1)-k*u(1,1))*m^-1;
        for i = 1:N %recorriendo particiones
            times(i+1,1) = dtj*i; %guardando el tiempo
            %predictores
            pu(i+1,1) = u(i,1) + dtj*v(i,1) + ((dtj^2)/2)*(1-2*beta)*ac(i,1); 
            pv(i+1,1) = v(i,1) + (1-gamma)*dtj*ac(i,1);
            %correctores
            u(i+1,1) = pu(i+1,1) + beta*(dtj^2)*ac(i+1,1);
            v(i+1,1) = pv(i+1,1) + gamma*dtj*ac(i+1,1);
            ac(i+1,1) = (p(dtj*i) - c*v(i+1,1) - k*u(i+1,1))/(m+gamma*dtj*c + beta*(dtj^2)*k);
            %comprobar dt_max < Omega_critico/wn_max
        end

        figure
        hold on
        plot(times,u);
        title('Desplazamientos');
        xlabel('tiempo [sec]');
        ylabel('Desplazamiento [in]');
        hold off
        
        figure
        hold on
        plot(times,v);
        title('Velocidades')
        xlabel('tiempo [sec]');
        ylabel('Velocidad [in/sec]');
        hold off
        
        figure
        hold on
        plot(times,ac);
        title('Aceleraciones')
        xlabel('tiempo [sec]');
        ylabel('Aceleración [in/sec^2]');
        hold off
    end
end
% //Tomando apuntes
% Estabilidad
% --Incondicionalmente estables
% =>  2*beta>= gamma >= 1/2
% --Condicionalmente estable
% =>  gamma >= 1/2   y   beta <= gamma/2      dt máximo < Omegacritco/wnmax
%                                              (dt es el paso)
%___________________________________________________________________________
% Métodos                      | beta       | gamma  |establidad
%______________________________|____________|________|____________________________
% Aceleración promedio         |1/4         |1/2     |incond.estable
%       (regla trapezoidal)    |            |        |
%______________________________|____________|________|______________________________
% Aceleración lineal           |            |        |
%       (explícito)            |1/6         |1/2     |cond Estable (Omega_critico = 2 *sqrt(3)
%______________________________|____________|________|____________________________
%
