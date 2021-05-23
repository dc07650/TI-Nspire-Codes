clear all; close all; clc;

A = [0 -1;
     1 1];
B = [0;
     1];
Q = [1 0;
     0 0];
R = [1];
lambda = eig(A);
disp(lambda);
for i = 1: length(lambda)
    if abs(lambda(i)) > 1
       disp("System unstable, redefine state matrix.")
       return
    end
end
disp("Initiating discrete lqr simulation");

P = 0;
while true
    oldP = P;
    K = -inv(R+B'*P*B)*B'*P*A;
    
    P = Q + K'*R*K + (A+B*K)'*P*(A+B*K);
     
    if abs(P - oldP) < 0.000001
        break;
    end
end
disp("P = ");
disp(P);
K = -inv(R+B'*P*B)*B'*P*A
disp("K = ");
disp(K);
disp("LQR: K = ");
disp(lqr(A, B, Q, R))