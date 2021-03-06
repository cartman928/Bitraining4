%% Initialize Parameters

clc
clear

alpha = 0;    %coefficient for block fading model
beta = 0.8^2;  %attenuation loss from non-direct antennas
n0 = 10^(-2);    %noise variance

iternums = 1:50; % number of iterations
N_Realizations = 2;

C1 = zeros(N_Realizations, length(iternums));
C2 = zeros(N_Realizations, length(iternums));
C3 = zeros(N_Realizations, length(iternums));

%% Start Loop
for Realization = 1 : N_Realizations
    Realization
        
    %% Random Channels
    H11 = (randn(2,2)+1i*randn(2,2))/sqrt(2);
    H22 = (randn(2,2)+1i*randn(2,2))/sqrt(2);
    H33 = (randn(2,2)+1i*randn(2,2))/sqrt(2);
    
    H12 = (randn(2,2)+1i*randn(2,2))/sqrt(2/beta); 
    H13 = (randn(2,2)+1i*randn(2,2))/sqrt(2/beta); 
    H21 = (randn(2,2)+1i*randn(2,2))/sqrt(2/beta);
    H23 = (randn(2,2)+1i*randn(2,2))/sqrt(2/beta); 
    H31 = (randn(2,2)+1i*randn(2,2))/sqrt(2/beta); 
    H32 = (randn(2,2)+1i*randn(2,2))/sqrt(2/beta);
    
    %Backward Channel
    Z11 = H11';
    Z22 = H22';
    Z33 = H33';
    
    Z12 = H21';
    Z13 = H31';
    Z21 = H12';
    Z23 = H32';
    Z31 = H13';
    Z32 = H23';   
    
    %% one iteration per block
    g1 = rand(2, 1) + 1i*rand(2, 1);    
    g2 = rand(2, 1) + 1i*rand(2, 1);
    g3 = rand(2, 1) + 1i*rand(2, 1);
    g1/norm(g1);
    g2/norm(g2);
    g3/norm(g3);
 
    v11 = zeros(2, 1); 
    v12 = zeros(2, 1);
    v13 = zeros(2, 1); 
    v21 = zeros(2, 1); 
    v22 = zeros(2, 1);
    v23 = zeros(2, 1); 
    v31 = zeros(2, 1); 
    v32 = zeros(2, 1); 
    v33 = zeros(2, 1); 
    
    for numiters = 1:length(iternums)

        %% bi-directional training
            
            %%Backward Training: sudo-LS Algorithm
            %abs(error)<2*10^(-4)
            [v11, v12, v13, v21, v22, v23, v31, v32, v33] = MSE_b(Z11, Z12, Z13, Z21, Z22, Z23, Z31, Z32, Z33, g1, g2, g3, n0);
            %[v11, v12, v13]
            %[v21, v22, v23]
            %[v31, v32, v33]
            %Power = [norm(v11)^2+norm(v12)^2+norm(v13)^2 norm(v21)^2+norm(v22)^2+norm(v23)^2 norm(v31)^2+norm(v32)^2+norm(v33)^2]
            %{
            v12 = [0;0];
            v13 = [0;0];
            v21 = [0;0];
            v23 = [0;0];
            v31 = [0;0];
            v32 = [0;0];
            v11 = v11/norm(v11);
            v22 = v22/norm(v22);
            v33 = v33/norm(v33);
            %}
          

            %%Forward Training: LS Algorithm
            [g1, g2, g3] = MSE_f(H11, H12, H13, H21, H22, H23, H31, H32, H33, v11, v12, v13, v21, v22, v23, v31, v32, v33, n0);
            %norm(g1)^2
            %norm(g2)^2

            SINR1 = norm(g1'*(H11*v11+H12*v21+H13*v31))^2/(norm(g1'*(H11*v12+H12*v22+H13*v32))^2+norm(g1'*(H11*v13+H12*v23+H13*v33))^2+n0*g1'*g1);
            SINR2 = norm(g2'*(H21*v12+H22*v22+H23*v32))^2/(norm(g2'*(H21*v11+H22*v21+H23*v31))^2+norm(g2'*(H21*v13+H22*v23+H23*v33))^2+n0*g2'*g2);
            SINR3 = norm(g3'*(H31*v13+H32*v23+H33*v33))^2/(norm(g3'*(H31*v12+H32*v22+H33*v32))^2+norm(g3'*(H31*v12+H32*v22+H33*v32))^2+n0*g3'*g3);
            C1(Realization, numiters) = abs(log2(1+SINR1));
            C2(Realization, numiters) = abs(log2(1+SINR2));
            C3(Realization, numiters) = abs(log2(1+SINR3));
            end           
    
end



%% Plot C(bits/channel)
%figure
%hold on

p1=plot(iternums, mean(C1)+mean(C2)+mean(C3),'--');

axis([1 numiters 0 20])

xlabel('Number of iterations')
ylabel('C(bits/channel)')
title('Simple Receivers')
%}
