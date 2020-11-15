% This Matlab script is the function "functionExpectations_ls" used to generate
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

function [B,C,D,Pl,Ik_coh, Ik_noncoh]=...
    functionExpectations_ls(HMean,betaa,K,L,N,rhoP,tau_p,sigma2)
%Number of UEs sharing the same pilot sequence
KK = K/tau_p;

%Prepare to store statistical matrices
R = zeros(N,N,L,K);
Psi = zeros(N,N,L,K);


for l = 1:L
    for k = 1:K
        R(:,:,l,k) = HMean(:,l,k)*HMean(:,l,k)'+betaa(l,k)*eye(N);
    end
end

for l = 1:L
    for pp = 1:tau_p
        Psiprep = rhoP*tau_p*sum(R(:,:,l,(pp-1)*KK+1:pp*KK),4)+sigma2*eye(N);
        for kk = 1:KK
            
            k = (pp-1)*KK+kk;
            Psi(:,:,l,k) = Psiprep;
        end
    end
end

%Prepare to store expectations in the paper
B = zeros(L,K);
C = zeros(L,L,K,K);
D = zeros(L,L,K);
Pl = zeros(L,K);
Ik_coh = zeros(K,L,L,K);
Ik_noncoh = zeros(K,L,K);

%Compute the expectations in the paper
for l = 1:L
    for k = 1:K
        B(l,k) = sqrt(tau_p*rhoP)*HMean(:,l,k)'*HMean(:,l,k)...
            +N*sqrt(tau_p*rhoP)*betaa(l,k);
        D(l,l,k) = sigma2*trace(Psi(:,:,l,k));
    end
    
end

for l = 1:L
    for k = 1:K
        for k2 = 1:K
            C(l,l,k,k2) = trace(Psi(:,:,l,k)*R(:,:,l,k2));
        end
        pp = ceil(k/KK);
        
        for kk = 1:KK
            
            k3 = (pp-1)*KK+kk;
            C(l,l,k,k3) = C(l,l,k,k3)+tau_p*rhoP*2*N*betaa(l,k3)...
                *real(HMean(:,l,k3)'*HMean(:,l,k3))+...
                tau_p*rhoP*N^2*betaa(l,k3)^2;
            
        end
    end
    for l2 = 1:L
        if l2~=l
            for k = 1:K
                pp = ceil(k/KK);
                for kk = 1:KK
                    k3 = (pp-1)*KK+kk;
                    
                    C(l,l2,k,k3) = tau_p*rhoP*...
                        real(HMean(:,l,k3)'*HMean(:,l,k3)+...
                        N*betaa(l,k3))*real(HMean(:,l2,k3)'*HMean(:,l2,k3)+...
                        N*betaa(l2,k3));
                end
            end
        end
    end
end

for l = 1:L
    for k = 1:K
        
        Pl(l,k) = trace(Psi(:,:,l,k));
        
    end
end

for k = 1:K
    rr = ceil(k/KK);
    
    for l = 1:L
        for kk = 1:KK
            for pp = 1:tau_p
                k2 = (pp-1)*KK+kk;
                
                Ik_coh(k,l,l,k2) = trace(Psi(:,:,l,k2)*R(:,:,l,k));
                Ik_noncoh(k,l,k2) = trace(Psi(:,:,l,k2)*R(:,:,l,k));
                
                if pp==rr
                    Ik_coh(k,l,l,k2) = Ik_coh(k,l,l,k2)+tau_p*rhoP*(2*N*betaa(l,k)*...
                        real(HMean(:,l,k)'*HMean(:,l,k))+N^2*betaa(l,k)^2);
                    Ik_noncoh(k,l,k2) = Ik_noncoh(k,l,k2)+tau_p*rhoP*(2*N*betaa(l,k)*...
                        real(HMean(:,l,k)'*HMean(:,l,k))+N^2*betaa(l,k)^2);
                    
                    for l2 = 1:L
                        if l2~=l
                            Ik_coh(k,l,l2,k2) = tau_p*rhoP*...
                                (real(HMean(:,l,k)'*HMean(:,l,k))+N*betaa(l,k))*...
                                (real(HMean(:,l2,k)'*HMean(:,l2,k))+N*betaa(l2,k));
                        end
                    end
                end
            end
        end
    end
end