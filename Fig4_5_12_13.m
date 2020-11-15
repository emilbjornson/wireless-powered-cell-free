% This Matlab script generates Figures 4, 5, 12, and 13 in the paper:
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

%Number of UEs
K = 20;

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

%Rectifier parameters for linear and non-linear energy harvesting model M1
a = 0.3929;
b = 0.01675;
c = 0.04401;

A_nonlin = (a*c-b)*1000;
B_nonlin = 1e+6*c;
C_nonlin = 1000*c^2;
effic = A_nonlin/C_nonlin;

%Prepare to save simulation results
%mmf: the proposed max-min fair power control
%pc: the benchmark power control scheme

%For saving SE values
LMMSE_mmf_coh_lin = zeros(K,nbrOfSetups);
LMMSE_mmf_noncoh_lin = zeros(K,nbrOfSetups);
LMMSE_mmf_coh_nonlin = zeros(K,nbrOfSetups);
LMMSE_mmf_noncoh_nonlin = zeros(K,nbrOfSetups);

LMMSE_pc_coh_lin = zeros(K,nbrOfSetups);
LMMSE_pc_noncoh_lin = zeros(K,nbrOfSetups);
LMMSE_pc_coh_nonlin = zeros(K,nbrOfSetups);
LMMSE_pc_noncoh_nonlin = zeros(K,nbrOfSetups);

LS_mmf_coh_lin = zeros(K,nbrOfSetups);
LS_mmf_noncoh_lin = zeros(K,nbrOfSetups);
LS_mmf_coh_nonlin = zeros(K,nbrOfSetups);
LS_mmf_noncoh_nonlin = zeros(K,nbrOfSetups);

LS_pc_coh_lin = zeros(K,nbrOfSetups);
LS_pc_noncoh_lin = zeros(K,nbrOfSetups);
LS_pc_coh_nonlin = zeros(K,nbrOfSetups);
LS_pc_noncoh_nonlin = zeros(K,nbrOfSetups);

%For saving downlink power coefficients 
DLMMSE_coh_lin = zeros(L,K,nbrOfSetups);
DLMMSE_noncoh_lin = zeros(L,K,nbrOfSetups);
DLMMSE_coh_nonlin = zeros(L,K,nbrOfSetups);
DLMMSE_noncoh_nonlin = zeros(L,K,nbrOfSetups);

DLS_coh_lin = zeros(L,K,nbrOfSetups);
DLS_noncoh_lin = zeros(L,K,nbrOfSetups);
DLS_coh_nonlin = zeros(L,K,nbrOfSetups);
DLS_noncoh_nonlin = zeros(L,K,nbrOfSetups);

%For saving uplink power coefficients 
ULMMSE_coh_lin = zeros(K,nbrOfSetups);
ULMMSE_noncoh_lin = zeros(K,nbrOfSetups);
ULMMSE_coh_nonlin = zeros(K,nbrOfSetups);
ULMMSE_noncoh_nonlin = zeros(K,nbrOfSetups);

ULS_coh_lin = zeros(K,nbrOfSetups);
ULS_noncoh_lin = zeros(K,nbrOfSetups);
ULS_coh_nonlin = zeros(K,nbrOfSetups);
ULS_noncoh_nonlin = zeros(K,nbrOfSetups);

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
    
    
    %Obtain the long-term statistical terms for SE and harvested energy
    %with LS-based channel estimation
    [B2,C2,D2,Pl2,Ik_coh2,Ik_noncoh2]=...
        functionExpectations_ls(HMean,channelGain_NLOS,K,L,N,rhoP,tau_p,sigma2);
    
    Pl2 = real(Pl2);
    Ik_coh2 = real(Ik_coh2);
    Ik_noncoh2 = real(Ik_noncoh2);
    
    %Obtain the SEs without prelog factor for the proposed MMF and 
    %the benchmark power control schemes, and power coefficients for the MMF 
    [rate_lmmse_noncoh_lin1,rate_lmmse_noncoh_lin2,...
        check_lmmse_noncoh_lin1,check_lmmse_noncoh_lin2, pow_lmmse_noncoh_lin, eta_lmmse_noncoh_lin]...
        =optimize_power_noncoh_lin(K,L,B1,C1,D1,Pl1,...
        Ik_noncoh1,rhoD,tau_u,tau_p,...
        tau_d,effic,rhoP);
    
    [rate_lmmse_noncoh_nonlin1,rate_lmmse_noncoh_nonlin2,...
        check_lmmse_noncoh_nonlin1,check_lmmse_noncoh_nonlin2, pow_lmmse_noncoh_nonlin, eta_lmmse_noncoh_nonlin]...
        =optimize_power_noncoh_nonlin(K,L,B1,C1,D1,Pl1,...
        Ik_noncoh1,rhoD,tau_u,tau_p,...
        tau_d,A_nonlin,B_nonlin,C_nonlin,rhoP);
    
    [rate_lmmse_coh_lin1,rate_lmmse_coh_lin2,...
        check_lmmse_coh_lin1,check_lmmse_coh_lin2, pow_lmmse_coh_lin, eta_lmmse_coh_lin]...
        =optimize_power_coh_lin(K,L,B1,C1,D1,Pl1,...
        Ik_coh1,rhoD,tau_u,tau_p,...
        tau_d,effic,rhoP);
    
    [rate_lmmse_coh_nonlin1,rate_lmmse_coh_nonlin2,...
        check_lmmse_coh_nonlin1,check_lmmse_coh_nonlin2, pow_lmmse_coh_nonlin, eta_lmmse_coh_nonlin]...
        =optimize_power_coh_nonlin(K,L,B1,C1,D1,Pl1,...
        Ik_coh1,rhoD,tau_u,tau_p,...
        tau_d,A_nonlin,B_nonlin,C_nonlin,rhoP);
    
    
    %%%%%%%%%%%%%%%%%%%%%
    
    
    [rate_ls_noncoh_lin1,rate_ls_noncoh_lin2,...
        check_ls_noncoh_lin1,check_ls_noncoh_lin2, pow_ls_noncoh_lin, eta_ls_noncoh_lin]...
        =optimize_power_noncoh_lin(K,L,B2,C2,D2,Pl2,...
        Ik_noncoh2,rhoD,tau_u,tau_p,...
        tau_d,effic,rhoP);
    
    [rate_ls_noncoh_nonlin1,rate_ls_noncoh_nonlin2,...
        check_ls_noncoh_nonlin1,check_ls_noncoh_nonlin2, pow_ls_noncoh_nonlin, eta_ls_noncoh_nonlin]...
        =optimize_power_noncoh_nonlin(K,L,B2,C2,D2,Pl2,...
        Ik_noncoh2,rhoD,tau_u,tau_p,...
        tau_d,A_nonlin,B_nonlin,C_nonlin,rhoP);
    
    [rate_ls_coh_lin1,rate_ls_coh_lin2,...
        check_ls_coh_lin1,check_ls_coh_lin2, pow_ls_coh_lin, eta_ls_coh_lin]...
        =optimize_power_coh_lin(K,L,B2,C2,D2,Pl2,...
        Ik_coh2,rhoD,tau_u,tau_p,...
        tau_d,effic,rhoP);
    
    [rate_ls_coh_nonlin1,rate_ls_coh_nonlin2,...
        check_ls_coh_nonlin1,check_ls_coh_nonlin2, pow_ls_coh_nonlin, eta_ls_coh_nonlin]...
        =optimize_power_coh_nonlin(K,L,B2,C2,D2,Pl2,...
        Ik_coh2,rhoD,tau_u,tau_p,...
        tau_d,A_nonlin,B_nonlin,C_nonlin,rhoP);
    
    %Collect all the binary variables showing the feasibility of the
    %problems
    check1=[check_lmmse_noncoh_lin1, ...
        check_lmmse_noncoh_nonlin1, ...
        check_lmmse_coh_lin1, ...
        check_lmmse_coh_nonlin1, ...
        check_ls_noncoh_lin1, ...
        check_ls_noncoh_nonlin1, ...
        check_ls_coh_lin1, ...
        check_ls_coh_nonlin1];
    
    check2=[check_lmmse_noncoh_lin2, ...
        check_lmmse_noncoh_nonlin2, ...
        check_lmmse_coh_lin2, ...
        check_lmmse_coh_nonlin2, ...
        check_ls_noncoh_lin2, ...
        check_ls_noncoh_nonlin2, ...
        check_ls_coh_lin2, ...
        check_ls_coh_nonlin2];
    
    
    %Store the results if all the problems are feasible
    if (min(check1)>0.5)&&(min(check2)>0.5)
        n = n+1;
        
        LMMSE_mmf_coh_lin(:,n) = prelogFactor*rate_lmmse_coh_lin1;
        LMMSE_mmf_noncoh_lin(:,n) = prelogFactor*rate_lmmse_noncoh_lin1;
        LMMSE_mmf_coh_nonlin(:,n) = prelogFactor*rate_lmmse_coh_nonlin1;
        LMMSE_mmf_noncoh_nonlin(:,n) = prelogFactor*rate_lmmse_noncoh_nonlin1;
        
        LMMSE_pc_coh_lin(:,n) = prelogFactor*rate_lmmse_coh_lin2;
        LMMSE_pc_noncoh_lin(:,n) = prelogFactor*rate_lmmse_noncoh_lin2;
        LMMSE_pc_coh_nonlin(:,n) = prelogFactor*rate_lmmse_coh_nonlin2;
        LMMSE_pc_noncoh_nonlin(:,n) = prelogFactor*rate_lmmse_noncoh_nonlin2;
        
        LS_mmf_coh_lin(:,n) = prelogFactor*rate_ls_coh_lin1;
        LS_mmf_noncoh_lin(:,n) = prelogFactor*rate_ls_noncoh_lin1;
        LS_mmf_coh_nonlin(:,n) = prelogFactor*rate_ls_coh_nonlin1;
        LS_mmf_noncoh_nonlin(:,n) = prelogFactor*rate_ls_noncoh_nonlin1;
        
        LS_pc_coh_lin(:,n) = prelogFactor*rate_ls_coh_lin2;
        LS_pc_noncoh_lin(:,n) = prelogFactor*rate_ls_noncoh_lin2;
        LS_pc_coh_nonlin(:,n) = prelogFactor*rate_ls_coh_nonlin2;
        LS_pc_noncoh_nonlin(:,n) = prelogFactor*rate_ls_noncoh_nonlin2;
        
        DLMMSE_coh_lin(:,:,n) = pow_lmmse_coh_lin;
        DLMMSE_noncoh_lin(:,:,n) = pow_lmmse_noncoh_lin;
        DLMMSE_coh_nonlin(:,:,n) = pow_lmmse_coh_nonlin;
        DLMMSE_noncoh_nonlin(:,:,n) = pow_lmmse_noncoh_nonlin;
        
        DLS_coh_lin(:,:,n) = pow_ls_coh_lin;
        DLS_noncoh_lin(:,:,n) = pow_ls_noncoh_lin;
        DLS_coh_nonlin(:,:,n) = pow_ls_coh_nonlin;
        DLS_noncoh_nonlin(:,:,n) = pow_ls_noncoh_nonlin;
        
        ULMMSE_coh_lin(:,n) = eta_lmmse_coh_lin;
        ULMMSE_noncoh_lin(:,n) = eta_lmmse_noncoh_lin;
        ULMMSE_coh_nonlin(:,n) = eta_lmmse_coh_nonlin;
        ULMMSE_noncoh_nonlin(:,n) = eta_lmmse_noncoh_nonlin;
        
        ULS_coh_lin(:,n) = eta_ls_coh_lin;
        ULS_noncoh_lin(:,n) = eta_ls_noncoh_lin;
        ULS_coh_nonlin(:,n) = eta_ls_coh_nonlin;
        ULS_noncoh_nonlin(:,n) = eta_ls_noncoh_nonlin;
    end
end

% Plot Fig. 4
figure
ecdf(vec(LMMSE_mmf_noncoh_nonlin))
hold on
ecdf(vec(LMMSE_mmf_noncoh_lin))
ecdf(vec(LMMSE_mmf_coh_nonlin))
ecdf(vec(LMMSE_mmf_coh_lin))
ecdf(vec(LMMSE_pc_noncoh_nonlin))
ecdf(vec(LMMSE_pc_noncoh_lin))
ecdf(vec(LMMSE_pc_coh_nonlin))
ecdf(vec(LMMSE_pc_coh_lin))

% Plot Fig. 5
figure
ecdf(min(LMMSE_mmf_noncoh_nonlin,[],1))
hold on
ecdf(min(LMMSE_mmf_noncoh_lin,[],1))
ecdf(min(LMMSE_mmf_coh_nonlin,[],1))
ecdf(min(LMMSE_mmf_coh_lin,[],1))
ecdf(min(LMMSE_pc_noncoh_nonlin,[],1))
ecdf(min(LMMSE_pc_noncoh_lin,[],1))
ecdf(min(LMMSE_pc_coh_nonlin,[],1))
ecdf(min(LMMSE_pc_coh_lin,[],1))
      
% Plot Fig. 12
figure
ecdf(vec(ULMMSE_noncoh_nonlin))
hold on
ecdf(vec(ULMMSE_coh_nonlin)) 
ecdf(vec(ULS_noncoh_nonlin)) 
ecdf(vec(ULS_coh_nonlin)) 
      
% Plot Fig. 13
figure
ecdf(vec(DLMMSE_noncoh_nonlin)) 
hold on
ecdf(vec(DLMMSE_coh_nonlin)) 
ecdf(vec(DLS_noncoh_nonlin)) 
ecdf(vec(DLS_coh_nonlin))     