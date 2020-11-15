% This Matlab script is the function "functionExampleSetup" used to generate
% the figures in the paper:
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

function [HMean,channelGain_NLOS,H] = functionExampleSetup(L,K,N,nbrOfRealizations)

%Set the length in meters of the total square area
squareLength = 100;

%Number of APs per dimension
nbrAPsPerDim = sqrt(L);

%Standard deviation of shadow fading in dB
sigma_sf_NLOS = 4; %for NLOS
sigma_sf_LOS = 3;  %for LOS

%Define the antenna spacing (in number of wavelengths)
antennaSpacing = 1/2; %Half wavelength distance

%Carrier frequency in GHz
fc = 3.4;

%Distance between APs in vertical/horizontal direction
interAPDistance = squareLength/nbrAPsPerDim;

%Deploy APs on the grid
locationsGridHorizontal = repmat(interAPDistance/2:interAPDistance:squareLength-interAPDistance/2,[nbrAPsPerDim 1]);
locationsGridVertical = locationsGridHorizontal';
APpositions = locationsGridHorizontal(:) + 1i*locationsGridVertical(:);


%Prepare to store fixed part of the normalized LOS channels, channel gain (in dB),
%distances and angles between APs and UEs, and
%probababilities of LoS 
HMeanNormalized = zeros(N,L,K);
channelGaindB = zeros(L,K);
distancesAP = zeros(L,K);
angleAP = zeros(L,K);
probLOS = zeros(L,K);

%Prepare to store the channel gains of LOS and NLOS parts, and the fixed
%part of the LOS channels, and phase-shifted LOS channels
channelGain_LOS = zeros(L,K);
channelGain_NLOS = zeros(L,K);
HMean = zeros(N,L,K);
HMeanPs = zeros(N,nbrOfRealizations,L,K);


%Rician factors \kappa
ricianFactor = db2pow(7+4*randn(L,K));

%Drop the UEs 
posX = rand(K,1)*squareLength;
posY = rand(K,1)*squareLength;
UEpositions = posX + 1i*posY;

%Compute the 3D distances by taking 4 meter vertical difference into account
%Compute the fixed parts of the normalized LOS channels
for k = 1:K
    distancesAP(:,k) = sqrt(4^2+(abs(APpositions-UEpositions(k,1))).^2);
    angleAP(:,k) = angle(UEpositions(k,1)-APpositions);
    for l = 1:L
        HMeanNormalized(:,l,k) = (exp(1i*2*pi.*(0:(N-1))*sin(angleAP(l,k))*antennaSpacing)); 
    end
end

for l = 1:L
    for k = 1:K
        if distancesAP(l,k)<=18
            probLOS(l,k) = 1;
        elseif distancesAP(l,k)<=37
            probLOS(l,k) = rand<exp(-(distancesAP(l,k)-18)/27);
        else
            probLOS(l,k) = rand<=0.5;
        end
    end
end

for l=1:L
    for k=1:K
        
        if probLOS(l,k)==1
            channelGaindB(l,k)= -16.9*log10(distancesAP(l,k))-32.8-20*log10(fc)...
                +sigma_sf_LOS*randn;
        else
            channelGaindB(l,k)= -43.3*log10(distancesAP(l,k))-11.5-20*log10(fc)...
                +sigma_sf_NLOS*randn;
        end
    end
    
end

%Go through all UEs and apply the channel gains to the mean vectors
for l = 1:L
    for k = 1:K
        
        if probLOS(l,k)==1 %The LoS Path exists, Rician Factor ~= 0
            channelGain_LOS(l,k) = ricianFactor(l,k)/(ricianFactor(l,k) +1 )*db2pow(channelGaindB(l,k));
            channelGain_NLOS(l,k) = 1/(ricianFactor(l,k) +1 )*db2pow(channelGaindB(l,k));
        else  %Pure NLoS case
            channelGain_LOS(l,k) = 0;
            channelGain_NLOS(l,k) = db2pow(channelGaindB(l,k));
        end
        %Scaling fixed part of the LOS channels 
        HMean(:,l,k) = sqrt(channelGain_LOS(l,k))*HMeanNormalized(:,l,k);
        
        %Apply random phase shifts to the fixed part of the LOS channels
        for tt = 1:nbrOfRealizations
            HMeanPs(:,tt,l,k) = exp(1i*2*pi*rand)*HMean(:,l,k);
        end
        
    end
end



%Generate the channel realizations
%Generate uncorrelated random variables
W = (randn(N,nbrOfRealizations,L,K)+1i*randn(N,nbrOfRealizations,L,K));

%Prepare to store channel realizations
H = zeros(N,nbrOfRealizations,L,K);
for l = 1:L
    for k = 1:K
        
        H(:,:,l,k) = sqrt(0.5*channelGain_NLOS(l,k))*W(:,:,l,k)+HMeanPs(:,:,l,k);
    end
end