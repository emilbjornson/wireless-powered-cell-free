% This Matlab script is the function "optimize_power_coh_nonlin" used to generate
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


function [rate1,rate2,check,check2,powers_record, etas_record]=...
    optimize_power_coh_nonlin(K,L,B,C,D,Pl,Ik,rhoD,tau_u,tau_p,...
    tau_d,Arec,Brec,Crec,rhoP)
%Prepare to store the results
rate1 = zeros(K,1); % SE values for the proposed MMF power control
rate2 = zeros(K,1); % SE values for the benchmark power control scheme
%Binary variables showing the feasibility of the methods
check = 0;
check2 = 0;
%Prepare to store the power coefficients
etas = zeros(K,1);
powers = zeros(L,K);
powers_record = zeros(L,K); % downlink power allocation coefficients
etas_record = zeros(K,1); % uplink power control coefficients

%The following scalings are only for preventing the numerical accuracy
%issues during CVX implementation and not change the SE or harvested energy expressions
scale1 = 1/norm(vec(B))*L*K;
B = B*scale1;
C = C*scale1*scale1;
D = D*scale1*scale1;

scale2 = 1/norm(vec(D))*K;
DD = D*scale2;
rhoP2 = rhoP*scale2;

scale3 = 1/norm(vec(Ik))*K/Brec;
Pl = Pl*scale3;
Ik = Ik*scale3;

%Modified bisection search algorithm
A = ones(L,K);
t_low = 0;
t_upp = 10;
t = (t_low+t_upp)/2;

diff = 100;
while  (diff>0.001)&&(t>0.001)
    
    zeta = zeros(K,K);
    zeta2 = zeros(K,1);
    zeta3 = zeros(K,1);
    for k = 1:K
        zeta2(k,1) = real(abs(A(:,k)'*B(:,k))^2);
        zeta3(k,1) = real(A(:,k)'*DD(:,:,k)*A(:,k));
        for k2 = 1:K
            zeta(k,k2) = real(A(:,k)'*C(:,:,k,k2)*A(:,k));
        end
    end
    cvx_begin sdp quiet
    variable eta(K,1)
    variable PP(L,L,K)
    variable eee(K,1)
    minimize sum(eta)
    subject to
    eta>=zeros(K,1);
    for k = 1:K
        (1+t)*eta(k,1)*zeta2(k,1)>=t*(zeta(k,:)*eta+zeta3(k,1));
        Brec/scale2*(tau_p*rhoP2+tau_u*eta(k,1))<=Brec*(tau_d*Arec/Brec)*(1-Crec*eee(k,1));
        norm([eee(k,1); Brec*vec(Ik(k,:,:,:))'*vec(PP)+Crec; sqrt(2) ])<= eee(k,1)+Brec*vec(Ik(k,:,:,:))'*vec(PP)+Crec;
        
    end
    for l = 1:L
        Pl(l,:)*vec(PP(l,l,:))/rhoD<=1;
    end
    for k = 1:K
        PP(:,:,k)==semidefinite(L)
    end
    cvx_end
    
    
    if cvx_status(1)=='S'
        maxx = 0;
        for l = 1:L
            maxx = max(maxx,Pl(l,:)*vec(PP(l,l,:)));
            
        end
        if maxx<rhoD
            PP = PP*rhoD/maxx;
        end
        powerr = zeros(L,K);
        PP2 = zeros(L,L,K);
        for l = 1:L
            powerr(l,:) = vec(PP(l,l,:));
            if Pl(l,:)*(powerr(l,:)).'>rhoD
                powerr(l,:) = powerr(l,:)*rhoD/(Pl(l,:)*(powerr(l,:)).');
            end
        end
        
        for k = 1:K
            PP2(:,:,k) = sqrt(powerr(:,k))*sqrt(powerr(:,k))';
        end
        for k = 1:K
            Pin = vec(Ik(k,:,:,:))'*vec(PP2);
            Eharv = Arec*tau_d*Pin/(Brec*Pin+Crec);
            if (tau_p*rhoP2+tau_u*eta(k,1))>scale2*Eharv
                eta(k,1) = (scale2*Eharv-tau_p*rhoP2)/tau_u;
            end
        end
        
        maxx = 0;
        for k = 1:K
            Pin = vec(Ik(k,:,:,:))'*vec(PP2);
            Eharv = Arec*tau_d*Pin/(Brec*Pin+Crec);
            maxx = max( eta(k,1)/((scale2*Eharv-tau_p*rhoP2)/tau_u),maxx);
        end
        if maxx<1
            eta = eta/maxx;
        end
        if min(eta)>=0
            
            minn = inf;
            
            for k = 1:K
                
                Mat1 = DD(:,:,k)-eta(k,1)*B(:,k)*B(:,k)';
                vec1 = sqrt(eta(k,1))*B(:,k);
                for k2 = 1:K
                    Mat1 = Mat1+eta(k2,1)*C(:,:,k,k2);
                end
                A(:,k) = inv(Mat1)*vec1;
                A(:,k) = A(:,k)/norm(A(:,k));
                
                minn = min(minn,real(abs(A(:,k)'*vec1)^2/(A(:,k)'*Mat1*A(:,k))));
            end
            etas = eta;
            powers = powerr;
            PP2s = PP2;
            t_low = minn;
            t_upp = 1.2*minn;
            
            
            
            check = 1;
        else
            t_upp = t;
        end
    else
        t_upp = t;
    end
    t = (t_low+t_upp)/2;
    diff = abs(t_low-t_upp)^2/abs(t)^2;
    
end

if check>0
    PP2 = PP2s;
    for l = 1:L
        if Pl(l,:)*(powers(l,:)).'>rhoD
            powers(l,:) = powers(l,:)*rhoD/(Pl(l,:)*(powers(l,:)).');
        end
    end
    
    
    for k = 1:K
        Pin = vec(Ik(k,:,:,:))'*vec(PP2);
        Eharv = Arec*tau_d*Pin/(Brec*Pin+Crec);
        if (tau_p*rhoP2+tau_u*etas(k,1))>scale2*Eharv
            etas(k,1) = (scale2*Eharv-tau_p*rhoP2)/tau_u;
        end
    end
    
    
    if min(etas)<0
        check = 0;
    end
end

if check>0
    
    zeta = zeros(K,K);
    for k = 1:K
        for k2 = 1:K
            zeta(k,k2) = real(A(:,k)'*C(:,:,k,k2)*A(:,k));
        end
    end
    for k = 1:K
        rate1(k,1) = log2(1+etas(k,1)*abs(A(:,k)'*B(:,k))^2/(zeta(k,:)*etas-etas(k,1)*abs(A(:,k)'*B(:,k))^2+real(A(:,k)'*DD(:,:,k)*A(:,k))));
    end
    powers_record = powers;
    etas_record = etas;
end

%The benchmark scheme
PP3 = zeros(L,L,K);
for l = 1:L
    powers(l,:) = 1./sqrt(Pl(l,:));
    powers(l,:) = powers(l,:)*rhoD/(Pl(l,:)*(powers(l,:)).');
end
for k = 1:K
    PP3(:,:,k) = sqrt(powers(:,k))*sqrt(powers(:,k))';
end
for k = 1:K
    Pin = vec(Ik(k,:,:,:))'*vec(PP3);
    Eharv = Arec*tau_d*Pin/(Brec*Pin+Crec);
    etas(k,1) = (scale2*Eharv-tau_p*rhoP2)/tau_u;
end
if min(etas)>=0
    check2 = 1;
    for k = 1:K
        
        Mat1 = DD(:,:,k)-etas(k,1)*B(:,k)*B(:,k)';
        vec1 = sqrt(etas(k,1))*B(:,k);
        for k2 = 1:K
            Mat1 = Mat1+etas(k2,1)*C(:,:,k,k2);
        end
        A(:,k) = inv(Mat1)*vec1;
        A(:,k) = A(:,k)/norm(A(:,k));
        
        rate2(k,1) = log2(1+real(abs(A(:,k)'*vec1)^2/(A(:,k)'*Mat1*A(:,k))));
    end
end