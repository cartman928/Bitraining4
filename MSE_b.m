function [v11, v12, v13, v21, v22, v23, v31, v32, v33] = MSE_b(Z11, Z12, Z13, Z21, Z22, Z23, Z31, Z32, Z33, g1, g2, g3, n0)
%update filters by MSE criterion

%MSE
R1 = (Z11*g1)*(Z11*g1)'+(Z12*g2)*(Z12*g2)'+(Z13*g3)*(Z13*g3)'+n0*eye(2);
R2 = (Z21*g1)*(Z21*g1)'+(Z22*g2)*(Z22*g2)'+(Z23*g3)*(Z23*g3)'+n0*eye(2);
R3 = (Z31*g1)*(Z31*g1)'+(Z32*g2)*(Z32*g2)'+(Z33*g3)*(Z33*g3)'+n0*eye(2);

P11 = Z11*g1;
P12 = Z12*g2;
P13 = Z13*g3;

P21 = Z21*g1;
P22 = Z22*g2;
P23 = Z23*g3;

P31 = Z31*g1;
P32 = Z32*g2;
P33 = Z33*g3;

v11 = R1\P11;
v12 = R1\P12;
v13 = R1\P13;

v21 = R2\P21;
v22 = R2\P22;
v23 = R2\P23;

v31 = R3\P31;
v32 = R3\P32;
v33 = R3\P33;

%Normalize

v11 = v11/norm();
v12 = R1\P12;
v13 = R1\P13;

v21 = R2\P21;
v22 = R2\P22;
v23 = R2\P23;

v31 = R3\P31;
v32 = R3\P32;
v33 = R3\P33;
end

    
    