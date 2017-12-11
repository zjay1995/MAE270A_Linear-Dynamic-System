function [gam,fre] = hinfnormc(ll,ul,tolerance,A,B,C,D)
I       = eye(length(D));
ww      = 0; %number of lo with good
q       = 0;
N       = 1;
Nmax    = 1000;
zer     = 1e-08;
check   = 0;
while N < Nmax
    gam     = (ll+ul)/2;
    Dl      = (gam)^2*I-D'*D;
    a       = A+B*inv(Dl)*D'*C;
    b       = -B*inv(Dl)*B';
    c       = C'*C+C'*D*inv(Dl)*D'*C;
    d       = -A'-C'*D*inv(Dl)*B';
    Aclp    = [a b; c d];
    reale   = real(eig(Aclp));
    image   = imag(eig(Aclp));
    eigd    = max(image);
    if(ul-ll)/2 < tolerance
        return
    end 
    
    for ii = (1: length(real(eig(Aclp))))
        if (abs(reale(ii)) < zer && image(ii) > zer) 
            check = 1;
        end 
    end
    N   = N+1;
    
    if check == 1
    ll  = gam;
    else
        ul = gam;
    end 
    check   = 0;
    term = (1+1j*eigd)/(1-1j*eigd);
    fre = log(term)/(ts*1j);
end 
