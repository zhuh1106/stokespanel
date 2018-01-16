%% mollifier function
function Psi = mypsi(dist)

Psi = NaN(size(dist));
H = @(r) exp(2*exp(-1./r)./(r-1));
Psi = H(dist);
Psi(dist>=1)=0;
Psi(dist<=0)=1;

end