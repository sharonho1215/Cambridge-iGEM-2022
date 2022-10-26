alpha = 10/60 
theta_1 = 0.05
gamma_p =  0.00038
eta = 0.018/60 
mu_2 = 1.5
mu_1 = 0.7
theta_2 = 0.005 
gamma_c = 0 
s = tf('s')
d = 0.5*theta_1



%disturbance no deg
gamma_c = 0
x_eq_rep = (mu_2 - mu_1)/(mu_1*theta_2)
z_1eq_rep = (alpha*mu_1*theta_2 - mu_2*gamma_p + gamma_p*mu_1)/(gamma_p*theta_1*(mu_2 - mu_1))
z_2eq_rep = gamma_p*mu_1*theta_1*(mu_2 - mu_1)/(eta*(alpha*mu_1*theta_2 - mu_2*gamma_p + gamma_p*mu_1))

A = [-gamma_p -theta_1*alpha/(theta_1*z_1eq_rep+1)^2 0;
    0 -eta*z_2eq_rep-gamma_c -eta*z_1eq_rep;
    -theta_2*mu_2/(theta_2*x_eq_rep+1)^2 -eta*z_2eq_rep -eta*z_1eq_rep-gamma_c]
B = [1;
    0;
    0]
C = [1 0 0]

D = 0

[b, a] = ss2tf(A, B, C, D)

G_no_deg_rep = tf(b, a)


opt = stepDataOptions('StepAmplitude',d);
[y,tOut] = step(G_no_deg_rep, opt)
total = y + x_eq_rep
plot(tOut,total)

%---------------------------------------------------------------------------------------------
%gamma_c = 0.00038
gamma_c = 0.00038
x_eq_rep_00038 = 209.88125683027403
z1_eq_rep_00038 = 21.794726966281953
z2_eq_rep_00038 = 105.79286828154224

A = [-gamma_p -theta_1*alpha/(theta_1*z1_eq_rep_00038+1)^2 0;
    0 -eta*z2_eq_rep_00038-gamma_c -eta*z1_eq_rep_00038;
    -theta_2*mu_2/(theta_2*x_eq_rep_00038+1)^2 -eta*z2_eq_rep_00038 -eta*z1_eq_rep_00038-gamma_c]
B = [1;
    0;
    0]
C = [1 0 0]

D = 0

[b, a] = ss2tf(A, B, C, D)

G_no_deg_rep_00038 = tf(b, a)


opt = stepDataOptions('StepAmplitude',d);
[y,tOut] = step(G_no_deg_rep_00038, opt)
total = y + x_eq_rep_00038
plot(tOut,total)

%---------------------------------------------------------------------------------------
%no deg integrator
q_1_rep = z_1eq_rep - z_2eq_rep
q_2_rep = z_2eq_rep + z_1eq_rep

A = [-gamma_p -theta_1/2*alpha/((theta_1/2*(q_1_rep+q_2_rep)+1)^2) -theta_1/2*alpha/((theta_1/2*(q_1_rep+q_2_rep)+1)^2);
    mu_2*theta_2/(theta_2*x_eq_rep+1)^2 -gamma_c 0;
    -theta_2*mu_2/(theta_2*x_eq_rep+1)^2 eta*q_1_rep -eta*q_2_rep-gamma_c]
B = [1;
    0;
    0]
C = [0 1 0]

D = 0

[b, a] = ss2tf(A, B, C, D)

integrator_rep_nodeg = tf(b, a)

opt = stepDataOptions('StepAmplitude',d)
[y,tOut] = step(integrator_rep_nodeg, opt)
total = y + q_1_rep
plot(tOut,total)

%---------------------------------------------------------------------------------------
%deg00038 integrator
gamma_c = 0.00038
x_eq_rep_00038 = 209.88125683027403
z1_eq_rep_00038 = 21.794726966281953
z2_eq_rep_00038 = 105.79286828154224

q_1_rep_00038 = z1_eq_rep_00038 - z2_eq_rep_00038
q_2_rep_00038 = z2_eq_rep_00038 + z1_eq_rep_00038

A = [-gamma_p -theta_1/2*alpha/((theta_1/2*(q_1_rep_00038+q_2_rep_00038)+1)^2) -theta_1/2*alpha/((theta_1/2*(q_1_rep_00038+q_2_rep_00038)+1)^2);
    mu_2*theta_2/(theta_2*x_eq_rep_00038+1)^2 -gamma_c 0;
    -theta_2*mu_2/(theta_2*x_eq_rep_00038+1)^2 eta*q_1_rep_00038 -eta*q_2_rep_00038-gamma_c]
B = [1;
    0;
    0]
C = [0 1 0]

D = 0

[b, a] = ss2tf(A, B, C, D)

integrator_rep_deg = tf(b, a)
opt = stepDataOptions('StepAmplitude',d)
[y,tOut] = step(integrator_rep_deg, opt)
total = y + q_1_rep_00038
plot(tOut,total)

%-----------------------------------------------------------------------------------------------------
%open loop - theta_2 no deg



%-----------------------------------------------------------------------------------------------------
%open loop - theta_2 deg 00038


%-----------------------------------------------------------------------------------------------------
%open loop - theta_1 no deg


%-----------------------------------------------------------------------------------------------------
