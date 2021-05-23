clear all; close all; clc;

A = [0 1;
     -8 -4];
B = [1;
     1];
Q = [1 0;
     0 0];
R = [0.000001];
disp("R >> 1: Then close to OL poles")
disp("R << 1: Then close to OL zeros")
disp(R);
lambda = eig(A);
disp(lambda);
for i = 1: length(lambda)
    if real(lambda(i)) > 0
       disp("System unstable, redefine state matrix.")
       return
    end
end
disp("Initiating discrete lqr simulation");

syms s;
disp([1 0]*inv(s*eye(length(lambda))-A)*B)
disp(eig(A-B*lqr(A, B, Q, R)))