% This Matlab script is the function "functionExpectations_lmmse_random_pilot_sequence"
%used to generate the figures in the paper:
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

function [B,C,D,Pl,Ik_coh, Ik_noncoh]=...
    functionExpectations_lmmse_random_pilot_sequence(HMean,betaa,K,L,N,eta,tau_p,sigma2,nbrOfRealizations,H)

% Random pilot generation
PilMatrix = randn(tau_p,K)+1i*randn(tau_p,K);
for k = 1:K
    PilMatrix(:,k) = PilMatrix(:,k)/norm(PilMatrix(:,k))*sqrt(tau_p);
end

PsiBar = kron(eye(N),PilMatrix);

SigmaX = zeros(N*K,N*K,L);

VectorChannel = zeros(N*K,nbrOfRealizations,L);
VectorEstimate = zeros(N*K,nbrOfRealizations,L);
for ll = 1:L
    
    for n = 1:N
        
        VectorChannel((n-1)*K+1:n*K,:,ll) = reshape(H(n,:,ll,:),[nbrOfRealizations,K]).';
        SigmaX((n-1)*K+1:n*K,(n-1)*K+1:n*K,ll) = diag(betaa(ll,:))+diag(abs(vec(HMean(n,ll,:))).^2);
        
        for m = [1:n-1 n+1:N]
            
            SigmaX((n-1)*K+1:n*K,(m-1)*K+1:m*K,ll) = diag(vec(HMean(n,ll,:).*conj(HMean(m,ll,:))));
            
        end
        
    end
    
end

Hhat = zeros(N,nbrOfRealizations,L,K);

for ll = 1:L
    
    ReceivedSignal = sqrt(eta)*PsiBar*VectorChannel(:,:,ll)...
        +sqrt(0.5*sigma2)*(randn(N*tau_p,nbrOfRealizations)+1i*randn(N*tau_p,nbrOfRealizations));
    
    VectorEstimate(:,:,ll) = sqrt(eta)*SigmaX(:,:,ll)*PsiBar'*...
        ((eta*PsiBar*SigmaX(:,:,ll)*PsiBar'+sigma2*eye(N*tau_p))\ReceivedSignal);
    
    for n = 1:N
        Hhat(n,:,ll,:) = (VectorEstimate((n-1)*K+1:n*K,:,ll)).';
    end
end


B = zeros(L,K);
C = zeros(L,L,K,K);
D = zeros(L,L,K);
Pl = zeros(L,K);
Ik_coh = zeros(K,L,L,K);
Ik_noncoh = zeros(K,L,K);

for l = 1:L
    for k = 1:K
        B(l,k) = trace(Hhat(:,:,l,k)'*H(:,:,l,k))/nbrOfRealizations;
        D(l,l,k) = sigma2*trace(Hhat(:,:,l,k)'*Hhat(:,:,l,k))/nbrOfRealizations;
        Pl(l,k) = trace(Hhat(:,:,l,k)'*Hhat(:,:,l,k))/nbrOfRealizations;
    end
        
end
for tt = 1:nbrOfRealizations
    
    
    for k = 1:K
        for k2 = 1:K
            Matr1 = reshape(Hhat(:,tt,:,k2),[N,L]);
            Matr2 = reshape(H(:,tt,:,k),[N,L]);

            TempVect = diag(Matr1'*Matr2);
            TempMat = (TempVect*TempVect')/nbrOfRealizations;
            C(:,:,k2,k) = C(:,:,k2,k)+TempMat;
            Ik_coh(k,:,:,k2) =  Ik_coh(k,:,:,k2)+reshape(TempMat,[1,L,L]);
                    
            Ik_noncoh(k,:,k2) = Ik_noncoh(k,:,k2)+reshape(diag(TempMat),[1,L]);
        end
    end
  
end
                     