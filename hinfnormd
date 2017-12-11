function [gam,fre] = hinfnormd(ll,ul,tolerance,A,B,C,D,ts)
Ia  = eye(length(A));
Ac  = -(Ia-A)/(Ia+A);
Bc  = sqrt(2)*inv(Ia+A)*B;
Cc  = sqrt(2)*C*inv(Ia+A);
Dc  = D-C*inv(Ia+A)*B;

[gam,fre] = hinfnormc(ll,ul,tolerance,Ac,Bc,Cc,Dc);
end
