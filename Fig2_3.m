% This Matlab script generates Figures 2 and 3 in the paper:
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



% Non-linear energy harvesting model M1 parameters
a1 = 0.3929;
b1 = 0.01675;
c1 = 0.04401;
    
A1 = (a1*c1-b1)*1000;
B1 = 1e+6*c1;
C1 = 1000*c1^2;


% Non-linear energy harvesting model M2 parameters
a2 = 2.463;    
b2 = 1.635;
c2 = 0.826;

A2 = (a2*c2-b2)*1000;
B2 = 1e+6*c2;
C2 = 1000*c2^2;

% Plot Fig. 2
x1 = linspace(0,0.0001,1000);
figure
% Non-linear energy harvesting model M1 with power values in mW
plot(x1*1000,A1*x1./(B1*x1+C1)*1000)
hold on
% Linear energy harvesting model with power values in mW
plot(x1*1000,A1*x1./C1*1000)


% Plot Fig. 3
x2 = linspace(0,0.001,1000);
figure
% Non-linear energy harvesting model M2 with power values in mW
plot(x2*1000,A2*x2./(B2*x2+C2)*1000)
