# functions needed to convert O2 from ml / l to kPa
# If units are in mg / l, first convert to ml / l using the 
# following conversion:
#o2 (ml / l) = o2 (mg / l) / 1.42905

library(gsw)
gsw_O2sol_SP_pt <- function(sal,pt) {
  x = sal
  pt68 = pt*1.00024
  y = log((298.15 - pt68)/(273.15 + pt68))
  
  a0 =  5.80871
  a1 =  3.20291
  a2 =  4.17887
  a3 =  5.10006
  a4 = -9.86643e-2
  a5 =  3.80369
  b0 = -7.01577e-3
  b1 = -7.70028e-3
  b2 = -1.13864e-2
  b3 = -9.51519e-3
  c0 = -2.75915e-7
  
  O2sol = exp(a0 + y*(a1 + y*(a2 + y*(a3 + y*(a4 + a5*y)))) + x*(b0 + y*(b1 + y*(b2 + b3*y)) + c0*x))
  return(O2sol)
}


calc_po2 <- function(depth, longitude, latitude, sal, temp, o2) {
  #O2 from trawl data is in ml/l - may need to be converted to umol/kg

  partial_molar_vol = 0.000032
  kelvin = 273.15
  boltz = 0.000086173324
  
  #calculate percent saturation for O2 - assumes  units of mL O2/L
  # Input:       S = Salinity (pss-78)
  #              T = Temp (deg C) ! use potential temp
  #depth is in meters
  #[umol/kg] = [ml/L]*44660/(sigmatheta(P=0,theta,S) + 1000)
  SA = gsw_SA_from_SP(sal,depth,longitude,latitude) #absolute salinity for pot T calc
  pt = gsw_pt_from_t(SA,temp,depth) #potential temp at a particular depth
  CT = gsw_CT_from_t(SA,temp,depth) #conservative temp
  sigma0 = gsw_sigma0(SA,CT)
  o2_umolkg = o2*44660/(sigma0+1000)
  
  
  O2_Sat0 = gsw_O2sol_SP_pt(sal,pt)
 
  #= o2satv2a(sal,pt) #uses practical salinity and potential temp - solubity at p =1 atm
  press = exp(depth*10000*partial_molar_vol/gas_const/(temp+kelvin))
  O2_satdepth = O2_Sat0*press
  
  #solubility at p=0
  sol0 = O2_Sat0/0.209
  sol_Dep = sol0*press
  po2 = o2_umolkg/sol_Dep
  return(po2)
}