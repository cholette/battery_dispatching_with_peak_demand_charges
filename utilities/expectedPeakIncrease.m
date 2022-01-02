function ecd = expectedPeakIncrease(k,uk,xk,recordedPeak,dk,sys)

C = sys.Ck(:,:,k);
D = sys.Dk(:,:,k);
E = sys.Ek(:,:,k);
H = sys.Hk(:,:,k);
Sw = sys.Swk(:,:,k);
Sv = sys.Svk(:,:,k);

sig = sqrt(H*Sw*H'+Sv);
mu = C*xk + E*dk + D*uk;
D2 = mu-recordedPeak;
alph = -D2/sig;
ecd = D2*(1-normcdf2(alph))+sig/sqrt(2*pi)*exp(-alph.^2/2);