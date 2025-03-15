% Modelo Dinamico Linealizado de Motor DC 
% Equipo 3 

% --- Parametros ---
    La    = 28.44e-3;  % Inductancia [H]
    Ra    = 1.27;      % Resistencia de armadura [Ohm]
    Kb    = 0.35;      % Constante de FEM [V·s/rad]
    Ka    = 0.35;      % Constante de torque [N·m/A]
    Jm    = 0.007;     % Inercia [kg·m^2]
    b     = 0.00173;   % Fricción viscosa [N·m·s/rad]
    tau_c = 0.0;       % Par de carga [N·m]

% --- Restricciones --- 
    vmin = -9.0; 
    vmax = 9.0; 

% --- Variables globales --- 
    global w_ref w_dot_ref w_dot;
    w_dot = 0; 
    w_ref = 20;         % Velocidad angular deseada 
    w_dot_ref = 0;      % Aceleracion angular 

% --- Variables de tiempo --- 
    t0 = 0; 
    tf = 60; 
    dt = 0.02; 
    tt = t0:dt:tf; 

% --- Condiciones iniciales --- 
    % x0 = [w, i_arm]
    x0 = [0, 0]; 

% --- Ganancias de PID --- 
    Kp = 5; 
    Ki = 0; 
    Kd = 3; 

% --- Arreglos para graficar --- 
    N = length(tt);
    u_control  = zeros(N,1);
    e_hist     = zeros(N,1);
    T_em_hist  = zeros(N,1);

[t, x] = ode45(@(t,state) motorDynamics(state, t, Kp, Kd, Ka, Jm, La, Kb, b, tau_c, Ra, vmin, vmax), tt, x0);

% Extraer soluciones
w_sol = x(:,1);
i_arm_sol = x(:,2);

for k = 1:N
    w_val = w_sol(k);
    i_arm_val = i_arm_sol(k);

    w_ref_val = w_ref_t(t(k));
    w_dot_ref_val = w_dot_ref_t(t(k));
    w_ddot_ref_val = w_ddot_ref_t(t(k));  % Aproximación numérica de la segunda derivada

    % Aceleración actual estimada
    w_dot_est = (1/Jm)*(Ka*i_arm_val - b*w_val - tau_c);

    % Error
    e = w_ref_val - w_val;
    e_dot = w_dot_ref_val - w_dot_est;

    % Control Auxiliar (u_aux): incorpora el término feedforward w_ddot_ref
    u_aux = w_ddot_ref_val + Kp*e + Kd*e_dot;
    
    % Entrada de control u(t) mediante linearización por retroalimentación:
    u_val = (Jm*La/Ka)*(u_aux - ((Ka/(Jm*La)*(-Ra*i_arm_val - Kb*w_val) - (b/(Jm^2))*(Ka*i_arm_val - b*w_val - tau_c))));
    u_val = max(vmin, min(vmax, u_val));

    % Guardar señales
    u_control(k) = u_val;
    e_hist(k)    = e;
    T_em_hist(k) = Ka * i_arm_val;  % Par electromagnético
end

%% --- Graficar en una cuadrícula 3x2 ---
figure;

% (1) Velocidad Angular + Setpoint
subplot(3,2,1);
plot(t, w_sol, 'LineWidth', 2); hold on;
plot(t, arrayfun(@w_ref_t, t), '--r', 'LineWidth', 2);
title('Velocidad Angular');
xlabel('Tiempo (s)');
ylabel('w (rad/s)');
legend('w','w_{ref}');
grid on;

% (2) Corriente de Armadura
subplot(3,2,2);
plot(t, i_arm_sol, 'LineWidth', 2);
title('Corriente de Armadura');
xlabel('Tiempo (s)');
ylabel('i_{arm} (A)');
grid on;

% (3) Entrada de Control (Voltaje)
subplot(3,2,3);
plot(t, u_control, 'LineWidth', 2);
title('Entrada de Control');
xlabel('Tiempo (s)');
ylabel('Voltaje (V)');
grid on;

% (4) Error
subplot(3,2,4);
plot(t, e_hist, 'LineWidth', 2);
title('Error');
xlabel('Tiempo (s)');
ylabel('e (rad/s)');
grid on;

% (5) Par Electromagnético
subplot(3,2,5);
plot(t, T_em_hist, 'LineWidth', 2);
title('Par Electromagnético');
xlabel('Tiempo (s)');
ylabel('T_{em} (N·m)');
grid on;

% (6) Subplot libre (opcional)
subplot(3,2,6);
text(0.3, 0.5, 'Nadota', 'FontSize', 20);
axis off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Funciones auxiliares

function w_val = w_ref_t(t)
    switch true
        case (t >= 0 && t < 5)
            w_val = 4 * t;
        case (t >= 5 && t < 10)
            w_val = 20;
        case (t >= 10 && t < 15)
            w_val = 20 + (t - 10);
        case (t >= 15 && t < 20)
            w_val = 25;
        case (t >= 20 && t < 25)
            w_val = 25 - 5 * (t - 20);
        case (t >= 25 && t < 30)
            w_val = 0;
        case (t >= 30 && t < 35)
            w_val = -5 * (t - 30);
        case (t >= 35 && t < 50)
            w_val = -25;
        case (t >= 50 && t < 55)
            w_val = -25 + 5 * (t - 50);
        case (t >= 55 && t <= 60)
            w_val = 0;
        otherwise
            w_val = NaN;  % Valor indefinido fuera del rango [0,60]
    end
end

function w_dot_val = w_dot_ref_t(t)
    switch true
        case (t >= 0 && t < 5)
            w_dot_val = 4;
        case (t >= 5 && t < 10)
            w_dot_val = 0;
        case (t >= 10 && t < 15)
            w_dot_val = 1;
        case (t >= 15 && t < 20)
            w_dot_val = 0;
        case (t >= 20 && t < 25)
            w_dot_val = -5;
        case (t >= 25 && t < 30)
            w_dot_val = 0;
        case (t >= 30 && t < 35)
            w_dot_val = -5;
        case (t >= 35 && t < 50)
            w_dot_val = 0;
        case (t >= 50 && t < 55)
            w_dot_val = 5;
        case (t >= 55 && t <= 60)
            w_dot_val = 0;
        otherwise
            w_dot_val = NaN;  % Valor indefinido fuera del rango [0,60]
    end
end

function w_ddot_val = w_ddot_ref_t(t)
    % Aproximación numérica de la segunda derivada mediante diferencia central
    dt_small = 1e-3; % Paso pequeño
    w_dot_plus = w_dot_ref_t(t + dt_small);
    w_dot_minus = w_dot_ref_t(t - dt_small);
    w_ddot_val = (w_dot_plus - w_dot_minus) / (2*dt_small);
end

function dstatedt = motorDynamics(state, t, Kp, Kd, Ka, Jm, La, Kb, b, tau_c, Ra, vmin, vmax)
    global w_dot; 
    w = state(1); 
    i_arm = state(2);  

    % Referencias 
    w_ref = w_ref_t(t); 
    w_dot_ref = w_dot_ref_t(t); 
    w_ddot_ref = w_ddot_ref_t(t); 

    % Calcular errores 
    e = w_ref - w; 
    e_dot = w_dot_ref - w_dot; 

    % Control Auxiliar (u_aux)
    u_aux = w_ddot_ref + Kp*e + Kd*e_dot; 
    
    % Entrada de control u(t) mediante linearización por retroalimentación:
    u = (Jm*La/Ka)*(u_aux - ((Ka/(Jm*La)*(-Ra*i_arm - Kb*w) - (b/(Jm^2))*(Ka*i_arm - b*w - tau_c))));
    u = max(vmin, min(vmax, u));

    % Dinámica del sistema
    di_arm = (1/La) * (u - Ra*i_arm - Kb*w); 
    w_dot = (1/Jm)*(Ka*i_arm - b*w - tau_c); 
    
    dstatedt = [w_dot ; di_arm]; 
end
