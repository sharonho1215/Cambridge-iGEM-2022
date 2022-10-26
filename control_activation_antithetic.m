%no deg

alpha = 10/60 
theta_1 = 0.005
gamma_p =  0.00038
eta = 0.018/60 
mu_2 = 1.5
mu_1 = 0.7
theta_2 = 0.0005 
gamma_c = 0 
s = tf('s')



x_eq = mu_1/theta_2
z1_eq = gamma_p * mu_1/(theta_1 * theta_2)
z2_eq = theta_1 * theta_2 / (eta * gamma_p)
d = 0.75*theta_1
d = 0.525

A = [-gamma_p theta_1 0;
    0 -eta*z2_eq-gamma_c -eta*z1_eq;
    theta_2 -eta*z2_eq -eta*z1_eq-gamma_c]
B = [1;
    0;
    0]
C = [1 0 0]

D = 0

[b, a] = ss2tf(A, B, C, D)

G_no_deg = tf(b, a)

G_no_deg_ss = ss(A, B, C, D)

opt = stepDataOptions('StepAmplitude',d);
[y,tOut] = step(G_no_deg, opt)
total = y + x_eq
plot(tOut,total)

%gamma_c = 0.00038
%------------------------------------------------------------------------------------------------------------

gamma_c = 0.00038
x_eq_00038 = 1339.114566237159
z1_eq_00038 = 101.772696976846
z2_eq_00038 = 21.66024249388398


A = [-gamma_p theta_1 0;
    0 -eta*z2_eq-gamma_c -eta*z1_eq;
    theta_2 -eta*z2_eq -eta*z1_eq-gamma_c]
B = [1;
    0;
    0]
C = [1 0 0]

D = 0

[b, a] = ss2tf(A, B, C, D)

G_00038 = tf(b, a)

opt = stepDataOptions('StepAmplitude',d);
[y,tOut] = step(G_00038, opt)
total = y + x_eq_00038
plot(tOut,total)

pole(G_00038)

%----------------------------------------------------------------------------
% no deg integrator
gamma_c = 0
q_1 = z1_eq - z2_eq
q_2 = z2_eq + z1_eq
A = [-gamma_p theta_1/2 theta_1/2;
    -theta_2 0 0;
    theta_2 eta*q_1 -eta*q_2]
B = [1;
    0;
    0]
C = [0 1 0]

D = 0

[b, a] = ss2tf(A, B, C, D)

integrator_no_deg = tf(b, a)

opt = stepDataOptions('StepAmplitude',d)
[y,tOut] = step(integrator_no_deg, opt)
total = y + q_1
plot(tOut,total)

%-----------------------------------------------------------------------------
%deg00038 integrator
q_1_00038 = z1_eq_00038 - z2_eq_00038
q_2_00038 = z2_eq_00038 + z1_eq_00038
A = [-gamma_p theta_1/2 theta_1/2;
    -theta_2 -gamma_c 0;
    theta_2 eta*q_1_00038 -eta*q_2_00038-gamma_c]
B = [1;
    0;
    0]
C = [0 1 0]

D = 0

[b, a] = ss2tf(A, B, C, D)

integrator_deg_00038 = tf(b, a)

opt = stepDataOptions('StepAmplitude',d)
[y,tOut] = step(integrator_deg_00038, opt)
total = y + q_1_00038
plot(tOut,total)



%-----------------------------------------------------------------------------------------------------------
%open loop - theta_2 no deg
gamma_c = 0
A = [-gamma_p theta_1, 0;
     0, -eta*z2_eq-gamma_c, -eta*z1_eq;
     0, -eta*z2_eq, -eta*z1_eq - gamma_c]
B = [0;
     0;
    -1]
C = [1, 0, 0]
D = 0

[b,a] = ss2tf(A, B, C, D)
open_theta2 = tf(b,a)
[re,im] = nyquist(open_theta2)

bodemag(1/(1+open_theta2)) %1.53dB
title('Sensitivity Function of \theta_2 (no degradation)')



unstable = open_theta2/(1+10.32133*open_theta2)
pole(unstable)
impulse(unstable)


%---------------------------------------------------------------------------------------------------------------
%open loop - theta 2 deg00038
gamma_c = 0.00038
A = [-gamma_p theta_1, 0;
     0, -eta*z2_eq_00038-gamma_c, -eta*z1_eq_00038;
     0, -eta*z2_eq_00038, -eta*z1_eq_00038 - gamma_c]
B = [0;
     0;
    -1]
C = [1, 0, 0]
D = 0

[b,a] = ss2tf(A, B, C, D)
open_theta2_deg = tf(b,a)
nyquist(open_theta2_deg, open_theta2)
legend('with degradation', 'without degradation')



unstable_nodeg = open_theta2/(1+0.0036*open_theta2)
pole(unstable_nodeg)
impulse(unstable_nodeg)
unstable_deg = open_theta2_deg/(1+0.00710955*open_theta2_deg)
pole(unstable_deg)
impulse(unstable_deg)

bodemag(1/(1+open_theta2_deg)) %1.52dB
title('Sensitivity Function of \theta_2 (degradation)')

bode(1/(1+open_theta2_deg), 1/(1+open_theta2))
%---------------------------------------------------------------------------------------------------------------
%open loop - theta_1 no deg
gamma_c = 0
A = [-gamma_p 0, 0;
     0, -eta*z2_eq-gamma_c, -eta*z1_eq;
     theta_2, -eta*z2_eq, -eta*z1_eq - gamma_c]
B = [-1;
     0;
    0]
C = [0, 1, 0]
D = 0

[b,a] = ss2tf(A, B, C, D)
open_theta1 = tf(b,a)
nyquist(open_theta1)

bodemag(1/(1+open_theta1))
title('Sensitivity Function of \theta_1 (no degradation)')
%----------------------------------------------------------------------------------------------------------------
%open loop - theta_1 deg

gamma_c = 0.00038
A = [-gamma_p 0, 0;
     0, -eta*z2_eq_00038-gamma_c, -eta*z1_eq_00038;
     theta_2, -eta*z2_eq_00038, -eta*z1_eq_00038 - gamma_c]
B = [-1;
     0;
    0]
C = [0, 1, 0]
D = 0

[b,a] = ss2tf(A, B, C, D)
open_theta1_deg = tf(b,a)

nyquist(open_theta1_deg, open_theta1)
bodemag(1/(1+open_theta1_deg))
title('Sensitivity Function of \theta_1 (degradation)')

bode(1/(1+open_theta1_deg), 1/(1+open_theta1))
