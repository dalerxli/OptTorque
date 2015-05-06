# express ZVAC as a "mathconst" type that expands to any desired precision
@Base.math_const ZVAC 376.73031346177065546819840042031930826862350835241 pi*BigFloat("119.9169832")
const pow = ^ # I got tired of converting pow(x,y) expressions into x^y
#/**********************************************************************/
#* Bessel                   
#**********************************************************************/
function getMN(X, L, aIn, Eps, Mu, Omega)
  #//Constants
  k = real(sqrt(Eps*Mu)*Omega); #//wavevector
  kz = cos(aIn*pi/180.0)*k;
  kt = sin(aIn*pi/180.0)*k;
  kpow2 = k^2;
  ktpow2 = kt^2;
  kzpow2 = kz^2;
  Z=ZVAC*sqrt(Mu/Eps) #relative wave impedance of exterior medium
  #////////////////////////////////////////
  # Coordinates (x,y,z) and (r, Phi, z)
    x = X[1]
    y = X[2]
    z = X[3]
    r=hypot(x,y)
    zpow2 = z^2
    rpow2 = (r^2)
    rpow3 = (r^3)
    rpow4 = (r^4)
  if(x==0.0 && y>=0.0)
    Phi=pi*0.5; 
  elseif(x==0.0 && y<0.0)
    Phi=-pi*0.5; 
  else
    Phi=atan2(y,x);
  end
  #////////////////////////////////////////	
  #//Bessel Functions Simplified
  #BesP =   jn(L,kt*r);
  #BesPm1 = jn(L-1,kt*r);
  #BesPm2 = jn(L-2,kt*r); 
  BesP =   besselj(L,kt*r);
  BesPm1 = besselj(L-1,kt*r);
  BesPm2 = besselj(L-2,kt*r); 
   
  #/// Set M, N accordingly 
  if(r!=0.0)
    #////////////////////////////////////////	
    #/// For Normal Cases 
    M0=exp(im*(L*Phi + kz*z))*((BesPm1*kt*y)/r - (BesP*L)/(im*x + y));
    M1=exp(im*(L*Phi + kz*z))*(-((BesPm1*kt*x)/r) + (BesP*L)/(x - im*y));
    M2=0.0;
    N0=(im*exp(im*(L*Phi + kz*z))*kz*((BesPm1*kt*x)/r - (BesP*L)/(x - im*y)))/k;
    N1=(exp(im*(L*Phi + kz*z))*kz*(-((BesP*L)/(x - im*y)) + (im*BesPm1*kt*y)/r))/k;
    N2=(BesP*exp(im*(L*Phi + kz*z))*ktpow2)/k;
  elseif(r==0.0&&L==1)
    #////////////////////////////////////////	
    #/// For L=1, r=0.0 
    M0=im*(0.5)*exp(im*kz*z)*kt;
    M1=-(exp(im*kz*z)*kt)/2.;
    M2=0.0;
    N0=(im*(0.5)*exp(im*kz*z)*kt*kz)/k;
    N1=-(exp(im*kz*z)*kt*kz)/(2.*k);
    N2=0.0;
  elseif(r==0.0&&L==0)
    #////////////////////////////////////////	
    #/// For L=0, r=0.0, only Nz is nonzero. 
    N2=(exp(im*kz*z)*ktpow2)/k;
  end
  #//end else following if(r==0)
  #/// For L>1, r=0.0, M=N=0. 
  return [M0,M1,M2,N0,N1,N2] # EH
end
