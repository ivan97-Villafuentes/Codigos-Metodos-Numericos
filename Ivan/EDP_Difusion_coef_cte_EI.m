% function EDP_Difusion_coef_cte_EI(N,M)
% Pb cont: u_t- k u_xx = 0, x\in [a,b], t\in [0,T]. 
% u(t,a)=0, u(t,b)=0, u(0,x)=U0.
% Esquema en espacio: Diferencias finitas centradas, soporte equiespaciado
% x de [a,b].
% Esquema en tiempo: Euler Implicito
% Dato de entrada: N= no. grados libertad (nodos interiores), % 
%                  M= no. etapas tiempo
% Datos salida:  plot en cada etapa de tiempo

function EDP_Difusion_coef_cte_EI(N,M) % N aprox. espacio, M aprox. tiempo

a=0;b=5; 
N1=N+1; h=(b-a)/N1 ; h2=h*h; % paso de espacio (uniforme)
X=(a+h):h:(b-h); X=X'; % nodos interiores en (a,b), vector columna de N componentes
X_amp=[a;X;b];  			% vector ampliado con los valores frontera

t=0; T=5; % tiempo inicial y final
dt=(T-t)/M; % paso de tiempo

% Datos estacionarios 
K=1.+0*X;  % conductividad

% montaje de matriz estacionaria (difusion)
Aleft_est=-K(1:N)*dt/h2;
Adiag=2*K(1:N)*dt/h2 + 1.;     
Aright_est=-K(1:N)*dt/h2;

% inicializacion en tiempo (Atencion compatibilidad con c.c.)
U0=X.*(5.-X);
U0_amp=[0;U0;0];

hold on;

plot(X_amp,U0_amp);


for k=1:M ;     % time loup 

A=spdiags([Aleft_est,Adiag,Aright_est],[-1 0 1],N,N);
% NOTA: con spdiags los vectores deben ser columnas con la misma dimension, y se
% pierden el primer número de la diagonal inferior y el último número de la diagonal superior 

B=U0; % segundo miembro, sin c.c.

U=A\B;  		% vector columna 

U_amp=[0;U;0]; 			% caso Dirichlet (vector columna)

pause(1)

plot(X_amp,U_amp) 		% dibuja la curva solucion aproximada

U0 = U;                     % update, vector columna
	
end       % end time loup
hold off;
end
