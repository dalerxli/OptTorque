# express ZVAC as a "mathconst" type that expands to any desired precision
@Base.math_const ZVAC 376.73031346177065546819840042031930826862350835241 pi*BigFloat("119.9169832")

const pow = ^ # I got tired of converting pow(x,y) expressions into x^y

#/**********************************************************************/
#* Gauss- Laguerre                                                    */
#**********************************************************************/
function getLG(X, P,L, w0, Eps,Mu,Omega)
    rt2=sqrt(2)	
    w0pow2 = (w0^2)
    w0pow3 = (w0^3)
    w0pow4 = (w0^4)
    
    # Coordinates (x,y,z) and (r, Phi, z)
    x = X[1]
    y = X[2]
    z = X[3]
    r=hypot(x,y)
    zpow2 = z^2
    rpow2 = (r^2)
    rpow3 = (r^3)
    rpow4 = (r^4)
    Phi=atan2(y,x)

    #/////////////////////////////////////////////////////////////////
    #Constants
    k = real(sqrt(Eps*Mu)*Omega) # wavevector
    kpow2 = k^2
    Z=ZVAC*sqrt(Mu/Eps) #relative wave impedance of exterior medium
    #/////////////////////////////////////////////////////////////////
    Labs = abs(L) # |L| in (double)
    Labspow2 = Labs^2
    Lpow2 = L^2
    #/////////////////////////////////////////////////////////////////
    if (L==0)
        Cnorm=sqrt(2.0/(pi*w0pow2))
    else
        Cnorm = sqrt(2.0*gamma(P+1.0)/(pi*w0pow2*gamma(P+Labs+1.0)))
    end
    #//////////////////////////////////////////////////////////////////
    zR =k*w0pow2/2.0 # Rayleigh Range
    zRpow2 = zR^2
    zRpow3 = zR^3
    #//////////////////////////////////////////////////////////////////
    w = w0*sqrt(zpow2/zRpow2+1.0) # beam waist radius w[z]
    wpow2 = (w^2)
    wpow3 = (w^3)
    wpow4 = (w^4)
    wpow5 = (w^5)
    wpow7 = (w^7)
    RR = z+zRpow2/z                 # Radius of curvature R[z]
    RRinv = z/(zpow2 + zRpow2)      # 1/R[z] 
    RRinvpow2 = RRinv^2
    zeta = atan2(z,zR)              # Gouy phase shift
    Rho = rt2*r/w                   # Rho[x,y,z]
    Rhopow2 = (Rho^2)
    Rhopow4 = (Rho^4)
    
    #////////////////////////////////////////////////////////////////
    #Laguerre Functions Simplified
    if (P==0)
        LagP=1.0
        LagPm1=0.0 
        LagPm2=0.0
    elseif (P==1)
        LagP= 1.0 + abs(L) - Rhopow2
        LagPm1=-1.0 
        LagPm2=0.0
    elseif(P==2)
        LagP=  lm_polynomial(   P,  abs(L), Rhopow2) 
        LagPm1=lm_polynomial(-1+P,1+abs(L), Rhopow2) 
        LagPm2=1.0 
    else
        LagP=  lm_polynomial(   P,  abs(L), Rhopow2) 
        LagPm1=lm_polynomial(-1+P,1+abs(L), Rhopow2) 
        LagPm2=lm_polynomial(-2+P,2+abs(L), Rhopow2) 
    end
    #////////////////////////////////////////////////////////////////
    #////////////////////////////////////////////////////////////////
    # cdouble uG0, M0L1[3], N0L1[3] # values at r=0 
  uG0=(exp(im*k*z - im*zeta)*w0)/w
  #////////////////////////////////////////////////////////////////
  #////////////////////////////////////////////////////////////////
  uG=(exp(im*(-0.5)*k*rpow2*RRinv - rpow2/wpow2 + im*k*z - im*zeta)*w0)/w
  #////////////////////////////////////////////////////////////////
  #////////////////////////////////////////////////////////////////
  if (r==0.0) 
      if (Labs==1.0)
          #/ Use M0L1, N0L1 for x=y=0, Abs[L]=1.  
          #/ if Abs[L]!=1, M=N=0 
          M0=(im*rt2*Cnorm*exp(im*(k*z - 2.0*(1.0 + P)*zeta))*(1.0 + P)*w0)/wpow2
          M1=-((rt2*Cnorm*exp(im*(k*z - 2.0*(1.0 + P)*zeta))*(1.0 + P)*w0)/wpow2)
          M2=0.0+0im
          N0=(im*rt2*Cnorm*exp(im*(k*z - 2.0*(1.0 + P)*zeta))*(1.0 + P)*w0)/wpow2
          N1=-((rt2*Cnorm*exp(im*(k*z - 2.0*(1.0 + P)*zeta))*(1.0 + P)*w0)/wpow2)
          N2=0.0+0im
      else
          #/ if Abs[L]!=1, M=N=0 
          M0=0.0+0im
          M1=0.0+0im
          M2=0.0+0im
          N0=0.0+0im
          N1=0.0+0im
          N2=0.0+0im
      end
  else #r!=0 
      M0=(Cnorm*exp(im*(L*Phi - (Labs + 2.0*P)*zeta))*pow(Rho,-1.0 + Labs)*uG*(-2.0*rt2*LagPm1*r*Rhopow2*w*y + LagP*(im*L*Rho*wpow2*x + r*(rt2*Labs*w + r*Rho*(-2.0 - im*k*RRinv*wpow2))*y)))/(rpow2*wpow2)
      M1=(Cnorm*exp(im*(L*Phi - (Labs + 2.0*P)*zeta))*pow(Rho,-1.0 + Labs)*uG*(2.0*rt2*LagPm1*r*Rhopow2*w*x + LagP*(-(rt2*Labs*r*w*x) + rpow2*Rho*(2.0 + im*k*RRinv*wpow2)*x + im*L*Rho*wpow2*y)))/(rpow2*wpow2)
      M2=0.0+0im
      if (z==0.0)
          N0=(im*(0.5)*Cnorm*exp(im*L*Phi - rpow2/w0pow2)*pow(Rho,-1.0 + Labs)*(2.0*(1.0 + Labs + 2.0*P)*(r*(2.0*LagP*r*Rho - rt2*Labs*LagP*w0 + 2.0*rt2*LagPm1*Rhopow2*w0)*x + im*LagP*L*Rho*w0pow2*y)*zR + im*k*LagP*L*Rho*w0pow2*y*(rpow2 - 2.0*zRpow2) + k*r*x*(-(rt2*Labs*LagP*w0*(rpow2 - 2.0*zRpow2)) + 2.0*rt2*LagPm1*Rhopow2*w0*(rpow2 - 2.0*zRpow2) + 2.0*LagP*r*Rho*(rpow2 - w0pow2 - 2.0*zRpow2))))/(k*rpow2*w0pow2*zRpow2)
          N1=(Cnorm*exp(im*L*Phi - rpow2/w0pow2)*pow(Rho,-1.0 + Labs)*(2.0*(1.0 + Labs + 2.0*P)*(LagP*L*Rho*w0pow2*x + im*r*(2.0*LagP*r*Rho - rt2*Labs*LagP*w0 + 2.0*rt2*LagPm1*Rhopow2*w0)*y)*zR + k*LagP*L*Rho*w0pow2*x*(rpow2 - 2.0*zRpow2) + im*k*r*y*(-(rt2*Labs*LagP*w0*(rpow2 - 2.0*zRpow2)) + 2.0*rt2*LagPm1*Rhopow2*w0*(rpow2 - 2.0*zRpow2) + 2.0*LagP*r*Rho*(rpow2 - w0pow2 - 2.0*zRpow2))))/(2.*k*rpow2*w0pow2*zRpow2)
          N2=(Cnorm*exp(im*L*Phi - rpow2/w0pow2)*pow(Rho,-2.0 + Labs)*(-4.0*LagP*rpow4*Rhopow2 + 4.0*rt2*rpow3*Rho*(Labs*LagP - 2.0*LagPm1*Rhopow2)*w0 + 2.0*rpow2*(-((-1.0 + Labs)*Labs*LagP) + 2.0*(LagP + LagPm1 + 2.0*Labs*LagPm1)*Rhopow2 - 4.0*LagPm2*Rhopow4)*w0pow2 + rt2*r*Rho*(-(Labs*LagP) + 2.0*LagPm1*Rhopow2)*w0pow3 + LagP*Lpow2*Rhopow2*w0pow4))/(k*rpow2*w0pow4)
      else
          N0=(-8.0*Cnorm*exp(im*(L*Phi - (Labs + 2.0*P)*zeta))*LagPm2*pow(Rho,2.0 + Labs)*uG*w0pow2*x*z)/(k*wpow4*zRpow2) + (Cnorm*exp(im*(L*Phi - (Labs + 2.0*P)*zeta))*LagPm1*pow(Rho,Labs)*((8.0*Labs*uG*w0pow2*x*zpow2)/wpow4 + (8.0*(1.0 + Labs)*uG*w0pow2*x*zpow2)/wpow4 + (4.0*rt2*Rho*uG*w0pow2*x*zpow2)/(r*wpow3) - (im*(4.0)*rt2*r*Rho*uG*(im*(-2.0) + k*RRinv*wpow2)*w0pow2*x*zpow2)/wpow5 - (im*(4.0)*rt2*L*Rho*uG*w0pow2*y*zpow2)/(r*wpow3) + (im*(4.0)*rt2*(Labs + 2.0*P)*Rho*RRinv*uG*x*zRpow3)/(r*w) + (2.0*rt2*Rho*uG*x*(2.0*wpow2*w0pow2*zpow2 + im*(2.0)*wpow4*zRpow2*(-(k*z) + RRinv*zR) + rpow2*(-4.0*w0pow2*zpow2 + im*k*RRinv*wpow4*(1.0 - 2.0*RRinv*z)*zRpow2)))/(r*wpow5)))/(2.*k*z*zRpow2) + (Cnorm*exp(im*(L*Phi - (Labs + 2.0*P)*zeta))*LagP*pow(Rho,Labs)*((-4.0*(-1.0 + Labs)*Labs*uG*w0pow2*x*zpow2)/(Rhopow2*wpow4) - (2.0*rt2*Labs*uG*w0pow2*x*zpow2)/(r*Rho*wpow3) + (2.0*rt2*Labs*r*uG*(2.0 + im*k*RRinv*wpow2)*w0pow2*x*zpow2)/(Rho*wpow5) + (im*(2.0)*rt2*Labs*L*uG*w0pow2*y*zpow2)/(r*Rho*wpow3) - (im*(2.0)*rt2*Labs*(Labs + 2.0*P)*RRinv*uG*x*zRpow3)/(r*Rho*w) - (2.0*(Labs + 2.0*P)*RRinv*uG*(im*(-2.0) + k*RRinv*wpow2)*x*zRpow3)/wpow2 - (2.0*L*(Labs + 2.0*P)*RRinv*uG*y*zRpow3)/rpow2 + (L*uG*y*(im*(2.0)*wpow2*w0pow2*zpow2 + 2.0*wpow4*zRpow2*(k*z - RRinv*zR) + rpow2*(im*(-4.0)*w0pow2*zpow2 + k*RRinv*wpow4*(-1.0 + 2.0*RRinv*z)*zRpow2)))/(rpow2*wpow4) + (rt2*Labs*uG*x*(rpow2*(4.0*w0pow2*zpow2 + im*k*RRinv*wpow4*(-1.0 + 2.0*RRinv*z)*zRpow2) + im*(2.0)*wpow2*(im*w0pow2*zpow2 + wpow2*zRpow2*(k*z - RRinv*zR))))/(r*Rho*wpow5) + (exp(-(rpow2/wpow2) - im*(0.5)*k*(rpow2*RRinv - 2.0*z) - im*zeta)*w0*x*(rpow2*(im*(-2.0) + k*RRinv*wpow2)*(im*(-4.0)*w0pow2*zpow2 + k*RRinv*wpow4*(-1.0 + 2.0*RRinv*z)*zRpow2) + 2.0*wpow2*((6.0 + im*k*RRinv*wpow2)*w0pow2*zpow2 + wpow2*zRpow2*(kpow2*RRinv*wpow2*z + im*k*(-(RRinv*wpow2) - 2.0*z + RRinvpow2*wpow2*(2.0*z + im*zR)) + im*(2.0)*RRinv*zR))))/wpow7))/(2.*k*z*zRpow2)
          N1=(-8.0*Cnorm*exp(im*(L*Phi - (Labs + 2.0*P)*zeta))*LagPm2*pow(Rho,2.0 + Labs)*uG*w0pow2*y*z)/(k*wpow4*zRpow2) + (Cnorm*exp(im*(L*Phi - (Labs + 2.0*P)*zeta))*LagPm1*pow(Rho,Labs)*((im*(4.0)*rt2*L*Rho*uG*w0pow2*x*zpow2)/(r*wpow3) + (8.0*Labs*uG*w0pow2*y*zpow2)/wpow4 + (8.0*(1.0 + Labs)*uG*w0pow2*y*zpow2)/wpow4 + (4.0*rt2*Rho*uG*w0pow2*y*zpow2)/(r*wpow3) - (im*(4.0)*rt2*r*Rho*uG*(im*(-2.0) + k*RRinv*wpow2)*w0pow2*y*zpow2)/wpow5 + (im*(4.0)*rt2*(Labs + 2.0*P)*Rho*RRinv*uG*y*zRpow3)/(r*w) + (2.0*rt2*Rho*uG*y*(2.0*wpow2*w0pow2*zpow2 + im*(2.0)*wpow4*zRpow2*(-(k*z) + RRinv*zR) + rpow2*(-4.0*w0pow2*zpow2 + im*k*RRinv*wpow4*(1.0 - 2.0*RRinv*z)*zRpow2)))/(r*wpow5)))/(2.*k*z*zRpow2) + (Cnorm*exp(im*(L*Phi - (Labs + 2.0*P)*zeta))*LagP*pow(Rho,Labs)*((im*(-2.0)*rt2*Labs*L*uG*w0pow2*x*zpow2)/(r*Rho*wpow3) - (4.0*(-1.0 + Labs)*Labs*uG*w0pow2*y*zpow2)/(Rhopow2*wpow4) - (2.0*rt2*Labs*uG*w0pow2*y*zpow2)/(r*Rho*wpow3) + (2.0*rt2*Labs*r*uG*(2.0 + im*k*RRinv*wpow2)*w0pow2*y*zpow2)/(Rho*wpow5) + (2.0*L*(Labs + 2.0*P)*RRinv*uG*x*zRpow3)/rpow2 - (im*(2.0)*rt2*Labs*(Labs + 2.0*P)*RRinv*uG*y*zRpow3)/(r*Rho*w) - (2.0*(Labs + 2.0*P)*RRinv*uG*(im*(-2.0) + k*RRinv*wpow2)*y*zRpow3)/wpow2 + (L*uG*x*(im*(4.0)*rpow2*w0pow2*zpow2 - im*(2.0)*wpow2*w0pow2*zpow2 + wpow4*zRpow2*(k*(-2.0*z + rpow2*(RRinv - 2.0*RRinvpow2*z)) + 2.0*RRinv*zR)))/(rpow2*wpow4) + (rt2*Labs*uG*y*(rpow2*(4.0*w0pow2*zpow2 + im*k*RRinv*wpow4*(-1.0 + 2.0*RRinv*z)*zRpow2) + im*(2.0)*wpow2*(im*w0pow2*zpow2 + wpow2*zRpow2*(k*z - RRinv*zR))))/(r*Rho*wpow5) + (exp(-(rpow2/wpow2) - im*(0.5)*k*(rpow2*RRinv - 2.0*z) - im*zeta)*w0*y*(rpow2*(im*(-2.0) + k*RRinv*wpow2)*(im*(-4.0)*w0pow2*zpow2 + k*RRinv*wpow4*(-1.0 + 2.0*RRinv*z)*zRpow2) + 2.0*wpow2*((6.0 + im*k*RRinv*wpow2)*w0pow2*zpow2 + wpow2*zRpow2*(kpow2*RRinv*wpow2*z + im*k*(-(RRinv*wpow2) - 2.0*z + RRinvpow2*wpow2*(2.0*z + im*zR)) + im*(2.0)*RRinv*zR))))/wpow7))/(2.*k*z*zRpow2)
          N2=(-8.0*Cnorm*exp(im*(L*Phi - (Labs + 2.0*P)*zeta))*LagPm2*pow(Rho,2.0 + Labs)*uG)/(k*wpow2) + (2.0*Cnorm*exp(im*(L*Phi - (Labs + 2.0*P)*zeta))*LagPm1*pow(Rho,Labs)*uG*(rt2*Rho*wpow2 + 2.0*r*(w + 2.0*Labs*w) + 2.0*rt2*rpow2*Rho*(-2.0 - im*k*RRinv*wpow2)))/(k*r*wpow3) + (Cnorm*exp(im*(L*Phi - (Labs + 2.0*P)*zeta))*LagP*pow(Rho,-2.0 + Labs)*uG*(-(rt2*Labs*r*Rho*wpow3) + Lpow2*Rhopow2*wpow4 + 2.0*rt2*Labs*rpow3*Rho*w*(2.0 + im*k*RRinv*wpow2) + rpow4*Rhopow2*pow(im*(-2.0) + k*RRinv*wpow2,2) + 2.0*rpow2*wpow2*(Labs - Labspow2 + Rhopow2*(2.0 + im*k*RRinv*wpow2))))/(k*rpow2*wpow4)
      end # else following if z == 0
  end # else following if(r==0.0)
          
  #/// For Azimuthal (TE) Polarization, let's use E=aM, H=aN. 
  #/// EH[0]=sum a[P,L]*M[P,L]
  #/// EH[1]=sum a[P,L]*N[P,L] 

  #/// For now, let's use a single mode of M for E and single mode of N for H. 
  return [M0,M1,M2,N0,N1,N2] # EH
end

#**********************************************************************/
#* Polynomials                            *****************************/
#**********************************************************************/

function lm_polynomial(n::Integer, m::Number, x::Number)
#****************************************************************************80
#    LM_POLYNOMIAL evaluates Laguerre polynomials Lm(n,m,x).
#  Recursion:
#    Lm(0,M,X)   = 1
#    Lm(1,M,X)   = (M+1-X)
#    if 2 <= N:
#      Lm(N,M,X)   = ( (M+2*N-1-X) * Lm(N-1,M,X)
#                   +   (1-M-N)    * Lm(N-2,M,X) ) / N
#  Licensing:
#    This code is distributed under the GNU LGPL license.
#  Modified:
#    10 March 2012
#  Author:
#    John Burkardt
#  Modified by Y. Eunnie Lee 2015, to input and output real numbers instead of arrays
#  Parameters:
#    Input, int N, the highest order polynomial to compute.
#    Input, int M, the parameter.  M must be nonnegative.
#    Input, double X, the evaluation points.
#    Output, double LM_POLYNOMIAL[MM*(N+1)], the function values.
    n < 0 && throw(ArgumentError("n = $n cannot be negative"))
    n == 0 && return one(x)
    vj2 = 1.0 #first term always 1.0 ; 
    vj1 = (m+1)-x
    for j = 2:n
        vj = (((m+2*j-1)-x)*vj1 + (-m-j+1)*vj2) / j
        vj2 = vj1
        vj1 = vj
    end
    return vj1
end
