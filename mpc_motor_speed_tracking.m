
clear  

clc
yalmip('clear');

Phorizon=10; % prediction horizon of MPC **variable

%% motor (plant) parameters initialization
R = 0.636;
Ld = 0.012;
Lq = 0.020;
ph = 0.088;
p = 5;
J = 1e-3;
B = 1.7e-3;
Ts = 7.8;
we0 = 1200;
id0 = 0.1;
iq0 = 8;
Tl = 5.4;


A = [-1*R/Ld, (Lq*we0)/Ld, (Lq*iq0)/Ld;
     -1*(Ld*we0)/Lq, -1*(R/Lq), -1*(Ld*id0)/Lq - ph/Lq;
     (3*(p.^2)*(Ld - Lq)*iq0)/(2*J), (3*(p.^2)*(Ld - Lq)*id0)/(2*J) + (3*(p.^2)*ph)/(2*J), -B/J];


B = [Ts/Ld, 0;
     0, Ts/Lq
     0,  0];

C = [1,  0,  0;
     0, 1, 0;
     0, 0, 1];

D = 0;

sys = ss(A,B,C,D);

Q = 100;
Q1=diag([1000;10;0.01]);
R = [1 0
     0 1];

R0 = 1200;
R1=[0.1;8;1200];
Ymin1 = [-500; -1; 0];
Ymax1 = [1; 400; 6000];
Umin = [-1000; 1000];
Umax = [-1000; 1000];

Ts=0.01;



%% control parameters initialization  
% c2d to discritize 
sysd=c2d(sys,Ts);

u_OPTIM  = sdpvar(2,1,Phorizon,'full'); % input variables at each time step in the prediction horizon 
U=sdpvar(2,1); % used to penalize the slew rate (smooth running of the actuators) 
d_u=sdpvar(2,1,Phorizon,'full'); % d_u, the change of the input variable from one time step to the other
x_OPTIM= sdpvar(3,1,Phorizon+1,'full'); % the state variable at each time stop along the prediction horizon

%% solver options: 
options = sdpsettings('solver','quadprog');
options = sdpsettings(options,'verbose',1 );
options = sdpsettings(options,'cachesolvers',1);
options;
options.quadprog; 


%% structure of the mpc problem at each time step
constraints=[];
objective = 0;

 for K = 1:Phorizon

   if K==1
       
       constraints=d_u(:,:,K)==u_OPTIM(:,:,K)-U;        
   else
       
    constraints=[constraints;d_u(:,:,K)==u_OPTIM(:,:,K)-(u_OPTIM(:,:,K-1))]; % constraints for penalizing slew rate
   end
      
   objective=objective...
+d_u(:,:,K)'*R*d_u(:,:,K)+(R1-x_OPTIM(:,:,K))'*Q1*(R1-x_OPTIM(:,:,K));
% objective has two terms: 1. d_u(:,:,K)'*R*d_u(:,:,K) to decrease changes in input (effort used to obtain tracking)
% 2. (C*x_OPTIM(:,:,K)-R0)'*Q*(C*x_OPTIM(:,:,K)-R0) to track the setpoint R0 
constraints=[constraints; x_OPTIM(:,:,K+1)==sysd.A*x_OPTIM(:,:,K)+...
         sysd.B*u_OPTIM(:,:,K); % constraints for system dynamics in discrete form

% URateMin<=d_u(:,:,K)<=URateMax; % constraints for the slew rate 
% Ymin1<=x_OPTIM(:,:,K)<=Ymax1; % constraints on output
% Umin <= u_OPTIM(:,:,K)<=Umax; % constraints on input
                   ];  
 end
 constraints=[constraints];
 objective=objective+(R1-x_OPTIM(:,:,K+1))'*Q1*(R1-x_OPTIM(:,:,K+1));
parameters_in = {x_OPTIM(:,:,1),U}; % initialized variables 
solutions_out = {u_OPTIM(:,:,:),x_OPTIM(:,:,2:end)}; % output from the mpc optimization problem
controller = optimizer(constraints,objective, options, parameters_in, solutions_out);  % to solver


%% prealocatting vectors and matrices

t_sim=100;
U_OPT_VECTOR    = zeros(2,Phorizon,length(t_sim));
X_OPT_VECTOR    = zeros(3,Phorizon,length(t_sim));
plot_speed=zeros(3,t_sim);
plot_u=zeros(2,t_sim);

    % loop constants
    counter = 0;
    count_diag = 0;       
    % structure 
   for t=1:t_sim
    

     if t==1 % initialize 
     x=[0;0;500];
      U=[100;100];
     else
         x=nxstate;
         U=u;
     end
 inputs = {x,U}; 
          [solutions,diagnostics] = controller{inputs}; 
 if diagnostics == 1 
  count_diag = count_diag + 1;
    if count_diag > 5
    error('The problem is infeasible');
    end
 end    
  
U_OPT_VECTOR(:,:,t) = solutions{1};
X_OPT_VECTOR(:,:,t) = solutions{2};

% Predictions obtained from solver:

u=U_OPT_VECTOR(:,1,t); % only the first optimal solution is used as input to the plant 
 

%% structure of control system below
%%     - 
%% R0 |------------>MPC--------->plant--------->|
%%    |+                                         | 
%%    |                                      nxstate
%%    |                                         |
%%    |                                         | 
%     |<-- - --- -------------------------------|                                     


nxstate=sysd.A*x+sysd.B*u;
plot_speed(:,t)=nxstate;
plot_u(:,t)=u;
 end

%% plots 

ref=1200*ones(1,t_sim);

figure(101)
subplot(3,1,1)
plot(plot_speed(1,:))
xlabel('Time steps (Ts=0.01)');
ylabel('Id (A)');
grid on

subplot(3,1,2)
hold on
plot(plot_speed(3,:),'b')
plot(ref,'r--')
legend('mpc','Ref')
xlabel('Time steps (Ts=0.01)');
ylabel('Speed (rpm)');
grid on

subplot(3,1,3)
plot(plot_speed(2,:))
xlabel('Time steps (Ts=0.01)');
ylabel('Iq (A)');
grid on

figure(102)
subplot(2,1,1)
plot(plot_u(1,:))
xlabel('Time steps (Ts=0.01)');
ylabel('Ud (V)');
grid on

subplot(2,1,2)
plot(plot_u(2,:))
xlabel('Time steps (Ts=0.01)');
ylabel('Uq (V)');
grid on

figure(103)
hold on
plot(plot_speed(3,:), 'b')
plot(ref, 'r--')
legend('mpc', 'Ref')
xlabel('Time steps (Ts=0.01)');
ylabel('Speed (rpm)');
grid on


%% RMSE

int=0;
for i=10:t_sim
    diff=(1200-plot_speed(i))^2;
    int=int+diff;
    RMSE=sqrt(int/(t_sim-10));
end

disp(RMSE)

%% simulation of open loop system

u = [0;0];

tspan = [0 0.4];
x0 = [1; 1; 1];
[t,x] = ode45(@(t,x) A*x + B*u, tspan, x0);
y = C*x.';
x1 = x(:,1);
x2 = x(:,2);
x3 = x(:,3);
figure;
plot(t, x(:,1), 'b-', t, x(:,2), 'r--', t, x(:,3)), 'g-';
xlabel('Time');
ylabel('States');
legend('I_d', 'I_q', 'W_e');
grid on 



    
    

