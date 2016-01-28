%% Initialize Parameters

clc
clear

alpha = 0;    %coefficient for block fading model
beta = 0.8^2;  %attenuation loss from non-direct antennas
w1 = 1;     %weight of MSE
w2 = 1;
w3 = 1;
n0 = 10^(-2);    %noise variance

iternums = 1:5; % number of iterations
N_Realizations = 100;

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
    
    %{
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
    %}
    
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
            for k1 = 1 : 1
            k1
            [v11, v12, v13, lambda1] = S_LS_User1_Brutal(H11, H12, H13, H21, H22, H23, H31, H32, H33, g1, g2, g3, v21, v22, v23, v31, v32, v33, n0, w1, w2, w3, 20000, 0.0001);
            [v21, v22, v23, lambda2] = S_LS_User2_Brutal(H11, H12, H13, H21, H22, H23, H31, H32, H33, g1, g2, g3, v11, v12, v13, v31, v32, v33, n0, w1, w2, w3, 20000, 0.0001);
            [v31, v32, v33, lambda3] = S_LS_User3_Brutal(H11, H12, H13, H21, H22, H23, H31, H32, H33, g1, g2, g3, v11, v12, v13, v21, v22, v23, n0, w1, w2, w3, 20000, 0.0001);
            %[v11, v12, v13]
            %[v21, v22, v23]
            %[v31, v32, v33]
            Power = [norm(v11)^2+norm(v12)^2+norm(v13)^2 norm(v21)^2+norm(v22)^2+norm(v23)^2 norm(v31)^2+norm(v32)^2+norm(v33)^2]
            Lambda = [lambda1 lambda2 lambda3]
            end


            %%Forward Training: LS Algorithm
            [g1, g2, g3] = LS(H11, H12, H13, H21, H22, H23, H31, H32, H33, v11, v12, v13, v21, v22, v23, v31, v32, v33, n0);
            %norm(g1)^2
            %norm(g2)^2

            SINR1 = norm(g1'*(H11*v11+H12*v21))^2/(norm(g1'*(H11*v12+H12*v22))^2+n0*g1'*g1);
            SINR2 = norm(g2'*(H21*v12+H22*v22))^2/(norm(g2'*(H21*v11+H22*v21))^2+n0*g2'*g2);
            C1(Realization, numiters, traininglength) = abs(log2(1+SINR1));
            C2(Realization, numiters, traininglength) = abs(log2(1+SINR2));
        
            
    end
            
    
end




%% Plot C(bits/channel)
figure
hold on

p1=plot(iternums, mean(C1(:,:,10))+mean(C2(:,:,10)),'--');

axis([1 numiters 0 12])
%}
