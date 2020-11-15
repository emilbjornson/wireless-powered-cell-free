% This Matlab script generates Figure 9 in the paper:
%
% Özlem Tugfe Demir and Emil Björnson,
% "Joint Power Control and LSFD for Wireless-Powered Cell-Free Massive MIMO,"
% IEEE Transactions on Wireless Communications, to appear.
%
% This is version 1.0 (Last edited: 2020-11-10)
%
% License: This code is licensed under the GPLv2 license. If you in any way
% use this code for research that results in publications, please cite our
% paper as described above.

%Number of APs 
L = 36;

%Number of UEs for two scenarios
KSet = [40 60];

%Number of antennas per AP 
N = 8;

%Select the number of setups with random UE locations
nbrOfSetups = 500;

%Pilot transmit power (W)
rhoP = 10^(-7);

%Total power limit per AP (W)
rhoD = 10/L;

%Compute noise power (in dBm)
noiseVariancedBm = -96;

%Compute noise power (dB)
noiseVariancedB = noiseVariancedBm-30;

%Compute noise power (W)
sigma2 = db2pow(noiseVariancedB);

%Select length of coherence block
tau_c = 200;

%Pilot Length
tau_p = 5;

%Number of downlink samples
tau_d = 25;

%Number of uplink samples
tau_u = tau_c-tau_d-tau_p;

%Prelog factor
prelogFactor = tau_u/tau_c;

%Rectifier parameters for non-linear energy harvesting model M1
a = 0.3929;
b = 0.01675;
c = 0.04401;

A_nonlin = (a*c-b)*1000;
B_nonlin = 1e+6*c;
C_nonlin = 1000*c^2;

%Prepare to save simulation results
%for saving SE values
LMMSE_mmf_coh_nonlin = zeros(max(KSet),nbrOfSetups,2);
LMMSE_mmf_noncoh_nonlin = zeros(max(KSet),nbrOfSetups,2);

%Go through each scenario
for scen = 1:2
    K = KSet(scen);
    %Initialize the number of setups
    n = 0;
    
    %% Go through all setups
    while n<nbrOfSetups
        
        %Obtain the fixed part of the LOS channels and channel gains of NLOS
        %parts
        [HMean,channelGain_NLOS,H] = functionExampleSetup(L,K,N,1);
        
        %Obtain the long-term statistical terms for SE and harvested energy
        %with LMMSE-based channel estimation
        [B1,C1,D1,Pl1,Ik_coh1,Ik_noncoh1]=...
            functionExpectations_lmmse(HMean,channelGain_NLOS,K,L,N,rhoP,tau_p,sigma2);
        
        Pl1 = real(Pl1);
        Ik_coh1 = real(Ik_coh1);
        Ik_noncoh1 = real(Ik_noncoh1);
        
        
        %Obtain the SEs without prelog factor for the proposed MMF 
      
        
        [rate_lmmse_noncoh_nonlin1,rate_lmmse_noncoh_nonlin2,...
            check_lmmse_noncoh_nonlin1,check_lmmse_noncoh_nonlin2, pow_lmmse_noncoh_nonlin, eta_lmmse_noncoh_nonlin]...
            =optimize_power_noncoh_nonlin(K,L,B1,C1,D1,Pl1,...
            Ik_noncoh1,rhoD,tau_u,tau_p,...
            tau_d,A_nonlin,B_nonlin,C_nonlin,rhoP);
        
        
        
        [rate_lmmse_coh_nonlin1,rate_lmmse_coh_nonlin2,...
            check_lmmse_coh_nonlin1,check_lmmse_coh_nonlin2, pow_lmmse_coh_nonlin, eta_lmmse_coh_nonlin]...
            =optimize_power_coh_nonlin(K,L,B1,C1,D1,Pl1,...
            Ik_coh1,rhoD,tau_u,tau_p,...
            tau_d,A_nonlin,B_nonlin,C_nonlin,rhoP);
        
       
        
        %Collect all the binary variables showing the feasibility of the
        %problems
        check1=[ check_lmmse_noncoh_nonlin1, check_lmmse_coh_nonlin1];
       
        
        
        %Store the results if all the problems are feasible
        if (min(check1)>0.5)
            n = n+1;
            
            LMMSE_mmf_coh_nonlin(1:K,n,scen) = prelogFactor*rate_lmmse_coh_nonlin1;
            LMMSE_mmf_noncoh_nonlin(1:K,n,scen) = prelogFactor*rate_lmmse_noncoh_nonlin1;
            
        end
    end
end

% Plot Fig. 9
figure
ecdf(vec(LMMSE_mmf_noncoh_nonlin(1:KSet(1),:,1)))
hold on
ecdf(vec(LMMSE_mmf_coh_nonlin(1:KSet(1),:,1)))
ecdf(vec(LMMSE_mmf_noncoh_nonlin(1:KSet(2),:,2)))
ecdf(vec(LMMSE_mmf_coh_nonlin(1:KSet(2),:,2)))