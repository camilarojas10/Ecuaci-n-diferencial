%Elaborado por: Andrea Camila Rojas Mondragón - 1913088
clear all;
close all;
clc;

%Se declaran los valores de h solicitados para las simulaciones
h1 = 0.5;
h2 = 0.25;
h3 = 0.125;

%Llamado de la función resultados con cada valor de h 
[h] = resultados(h1); 
[h] = resultados(h2); 
[h] = resultados(h3); 

% --------------------------------------------------------------------------

%Función que almacena todos los resultados y gráfica 
function [h] = resultados(h)
f = @(t, y) (1 + 4.* t).* sqrt(y); % Ecuación diferencial a resolver 
ft=@(t) (t/2+t.^2+1).^2; % Solución analítica de la ecuación diferencial 
yO = 1;  % Condición inicial
tspan = [0 3]; % Rango de simulación en 0 y 3
tabla = 1; % Variable utilizada para imprimir la tabla de resultados
t1 = tspan(1):(tspan(2) - tspan(1)) / 1000:tspan(2); %Vector "tiempo" con 1000 puntos
y_real= ft(t1);

fprintf('\n -> El valor de h para esta simulación es: %.3f \n', h);

[te, ye] = Eulode(f,tspan,yO,h,tabla); %Llamado de la función del método de Euler 
[emc_e]=cemc(te,ye); %Llamado de la función para calcular el error  

[t, y] = puntomedio(f, tspan, yO, h); %Llamado de la función del método de MidPoint
[emc_pm]=cemc(t,y);%Llamado de la función para calcular el error  

[tRK,yRK]= RK4M(f,tspan,yO,h);%Llamado de la función del método de Runge Kutta
[emc_rk]=cemc(tRK,yRK);%Llamado de la función para calcular el error  

%Concatenación de los valores del error para mostrarlos
eemc=[emc_e,emc_pm, emc_rk];
fprintf('\n Los errores medios cuadráticos en esta simulación son:\n ');
disp('    Euler          MidPoint         Runge Kutta        ');
fprintf(' %15.10f %15.10f %15.10f\n',eemc);

% Gráfica
figure;
plot(t1, y_real,'LineWidth', 4); % Gráfica de la solución analítica 
hold on 
plot(te, ye, 'ro', te,ye,'r','LineWidth', 2) %Gráfica de Euler.
hold on
plot(t, y, 'gpentagram', t, y, 'g','LineWidth', 2) %Gráfica de MidPoint.
hold on
plot(tRK, yRK, 'msquare', tRK, yRK, 'm', 'LineWidth', 2) %Gráfica de Runge Kutta.
hold on

%Simbología de las gráficas
xlabel('Tiempo'); ylabel('y(t)');
title('y(t)=(1 + 4t)*sqrt(y)(azul) vs Euler (rojo), MidPoint (verde) y Runge Kutta (magenta)');
hold off
end 

% --------------------------------------------------------------------------

%Función que calcula el error medio cuadrático 
function [emc] = cemc(a, b)
f = @(t, y) (1 + 4.* t).* sqrt(y);
ft=@(t) (t/2+t.^2+1).^2;
tspan = [0 3]; 
t1 = tspan(1):(tspan(2) - tspan(1)) / 1000:tspan(2); 
y_real= ft(t1); 
f_real=ft(a); %Se asignan los valores de la solución analítica  
n=length(a); % Se obtiene el tamaño para aplicarlo en la formula 
emc=(1/n)*sum((b-f_real).^2); %Se sustituyen los valores calculados  
end 

% -------------------------------------------------------------------------

%Función del métod de Euler (tomada de los códigos dados en clase)
function[t,y] = Eulode(f,tspan,yO,h,tabla)
    a= tspan(1);  
    b=tspan(2);   
    t=a:h:b; 
    y=zeros(1,length(t)); 

    y(1)=yO; 
    t(1)=a;  

    for i = 1:(length(y)-1) 
      y(i+1) =y(i)+h*f(t(i),y(i)); 
    end
    
    if tabla==1        
        fprintf('\n Euler \n');
        disp('Step Size      t                    y(t)');
        k = 1:length(t); 
        out= [k ; t; y] ;
        fprintf('%5d %15.10f       %15.10f\n',out)
    end
    fprintf('\n')
end

% --------------------------------------------------------------------

%Función del métod de MidPoint (tomada de los códigos dados en clase)
function [t, y] = puntomedio(f, tspan, yO, h)
    a = tspan(1);
    b = tspan(2);
    n = (b - a) / h;
    t = (a + h:h:b);
    k1 = feval(f, a, yO);
    k2 = feval(f, a + h/2, yO + k1/2 * h);
    y(1) = yO + k2 * h;
    for i = 1:n-1
        k1 = feval(f, t(i), y(i));
        k2 = feval(f, t(i) + h/2, y(i) + k1/2 * h);
        y(i + 1) = y(i) + k2 * h;
    end
    t = [a t];
    y = [yO y];
    fprintf('\n MidPoint \n');
    disp('   step         t                  y(t)')
    k = 1:length(t);
    out = [k; t; y];
    fprintf('%5d %15.10f       %15.10f\n', out);
    fprintf('\n')
end

% --------------------------------------------------------------------

%Función del métod de Runge Kutta de 4 orden (tomada de los códigos dados en clase)
function[t,y] = RK4M(f,tspan,yO,h)

    a= tspan(1); 
    b=tspan(2); 
    n=(b-a)/h ;
    t=a+h:h:b;
    
    k1=f(a,yO);
    k2=f(a+1/2*h,yO+1/2*k1*h);
    k3=f(a+1/2*h,yO+1/2*k2*h);
    k4=f(a+h,yO+k3*h);


    y(1) =yO+1/6*(k1+2*k2+2*k3+k4)*h;
    
    for i = 1:(length(t)-1)
        
        k1=f(t(i),y(i));
        k2=f(t(i)+1/2*h,y(i)+1/2*k1*h);
        k3=f(t(i)+1/2*h,y(i)+1/2*k2*h);
        k4=f(t(i)+h,y(i)+k3*h);

        y(i+1) =y(i)+1/6*(k1+2*k2+2*k3+k4)*h;
    
    end
    t = [a t]; 
    y = [yO y];
    fprintf('\nRunge Kutta de 4 orden\n')
    disp('   step         t              y(t)')
    k = 1:length(t); 
    out= [k ; t; y] ;
    fprintf('%5d %15.10f       %15.10f\n',out)
end

