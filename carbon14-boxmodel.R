# Box Diffusion carbon-cycle model for 14C only.
# This is a variant of the original Siegenthaler & Oeschger (1985) model that assumes a steady state carbon cycle
# No fractionation 
# State variables are the 14C/C ratios in all reservoirs
# These can be interpreted as deviation from a background steady state

# required library for solving differential equations
require(deSolve)

# time sequence 
time <- c(seq(0,100, by=1),seq(100,1000,by = 10),seq(1100,10000,by=100),seq(11000,200000,by=1000))

# conversion factors
gm = 1.2e-14 # gamma_mol; conversion factor from mol to GtC or PgC
gp = 2.123 # gamma_ppm; conversion factor from ppm to GtC or PgC

# atmosphere quantities
Pa0 = 270 # intial CO2 in the atmosphere (ppm) in 1751

# ocean parameters
K = 4000 # eddy diffusion constant in the deep sea (m2/y)
h1 = 120 # Delta_h thermocline (there are 10 intervals = 1200m)
h2 = 320.8333 # Delta_h deep ocean (there are 6 intervals = 2525m)
hm = 75  # average depth of mixed layer
Aoc = 3.62e+14 # ocean surface area
kex = 0.0522 # gas exchange coefficient (mol/m2/yr/muatm) (best estimate of Naegler, 2009)
Cm0 = 2 # CO2 concentration in the mixed layer (assumed to be constant with depth) (moles/m3)
akam = ( kex * Aoc * gm / gp)  # calculated from values defined above
tauam = 1/akam
akma = (kex * Pa0/ (hm * Cm0))

ak = K
akmd = ak / (hm * h1 * 2)
ak1 = ak / (h1 * h1)
ak2 = ak / (h2 * h2)
akv = (2 * ak) / (h1 * (h1 + h2))
akj = (2 * ak) / (h1 * h2)
akn = (2 * ak) / (h2 * (h1 + h2))

# land biosphere parameters
fab0 = 60 # net primary productivity (PgC/yr)
tauab = 60 # biosphere carbon turnover time (yr)
# tauam = 10 # 9-11 years: atmospheric turnovertime relative to mixed layer
akab = fab0/(gp * Pa0)
tauab = 1/akab
akba = 1/tauab

# conversion of Delta-14C to mol 14C: 14N = rs14 (1+Delta-14C)
# conversion of mol 14C to Delta-14C: Delta-14C = 14N/rs14 -1
# Delta-14C in absolute values, not in permil!

rs14 = 1.176e-12 #standard ratio of 14C/C (mol/mol, value of Naegler 2009)
alam14 = 1/8267 # 14C decay constant
Avo = 6.02e23 # Avogadro constant (at/mol)

# model parameters: a named vector
parameters <- c(
                hm = hm,
                gm = gm,
                gp = gp,
                gm = gm,
                Pa0 = Pa0,
                akab = akab,
                akba = akba,
                akam = akam,
                akma = akma,
                akmd = akmd, 
                ak1 = ak1,
                ak2 = ak2,
                akv = akv,
                akj = akj,
                akn = akn,
                rs14 = rs14,
                alam14 = alam14,
                Avo = Avo
)

# initial condition: a named vector
state <- c(Ra_14C = 0, Rb_14C = 0, Rm_14C = 0, 
           R1_14C = 0, R2_14C = 0, R3_14C = 0, R4_14C = 0, R5_14C = 0, 
           R6_14C = 0, R7_14C = 0, R8_14C = 0, R9_14C = 0, R10_14C = 0, R11_14C = 0, R12_14C = 0, R13_14C = 0,
           R14_14C = 0, R15_14C = 0, R16_14C = 0)

# One rectangular function (14C production in 14C at/yr) to force the model (i.e. 'perturbation' mode)
# F <- function(t) {
  # ifelse(t > 10 & t < 20, 10e26, 0)
# }

# Constant production for achieving steady state 
# (requires >100'000 year simulation because of decay constant of 14C)
# this value results in 1000 permil in the atmosphere
F <- function(t) {1.80841e26} # C14-atoms/yr

# Here you can prescribe transient changes in the gas exchange coefficient for perturbation experiments
# kex <- ifelse(t > 40 & t < 60, 0.64, 0.064)

# Function to calculate the value of the derivatives at each time value
# Use the names of the variables as defined in the vectors above
diffusion <- function(t, state, parameters){
  with(as.list(c(state, parameters)), {
    
    # Atmospheric production of 14C (in 14C atoms/yr); convert into change in atmospheric 14C/C ratio
    Q = F(t) * gm/(Avo * Pa0 * gp)
    
    
    #####################################
    ########## 14C equation set #########
    ######################################
    
    ####### atmosphere ########
    dRa_14C = Q + 
              akab * (Rb_14C - Ra_14C) +
              akam * (Rm_14C - Ra_14C) - 
              alam14 * Ra_14C
    
    ####### Land biosphere ######
    dRb_14C = akba * (Ra_14C - Rb_14C) -
              alam14 * Rb_14C
    
    ####### mixed layer ########
    dRm_14C = akma * (Ra_14C - Rm_14C)  + 
      akmd * ((-3*Rm_14C) + (4*R1_14C) - R2_14C) -
      alam14 * Rm_14C
    
    ####### thermocline ########
    dR1_14C = ak1 * (Rm_14C - (2* R1_14C) + R2_14C) - alam14 * R1_14C
    dR2_14C = ak1 * (R1_14C - (2* R2_14C) + R3_14C) - alam14 * R2_14C
    dR3_14C = ak1 * (R2_14C - (2* R3_14C) + R4_14C) - alam14 * R3_14C
    dR4_14C = ak1 * (R3_14C - (2* R4_14C) + R5_14C) - alam14 * R4_14C
    dR5_14C = ak1 * (R4_14C - (2* R5_14C) + R6_14C) - alam14 * R5_14C
    dR6_14C = ak1 * (R5_14C - (2* R6_14C) + R7_14C) - alam14 * R6_14C
    dR7_14C = ak1 * (R6_14C - (2* R7_14C) + R8_14C) - alam14 * R7_14C
    dR8_14C = ak1 * (R7_14C - (2* R8_14C) + R9_14C) - alam14 * R8_14C
    dR9_14C = ak1 * (R8_14C - (2* R9_14C) + R10_14C) - alam14 * R9_14C
    
    ####### Deep Ocean assuming a uniform grid #######
    dR10_14C = (akv * R9_14C) - (akj * R10_14C) + (akn * R11_14C) - alam14 * R10_14C 
    
    dR11_14C = ak2 * (R10_14C - (2* R11_14C) + R12_14C) - alam14 * R11_14C # deep ocean layer
    dR12_14C = ak2 * (R11_14C - (2* R12_14C) + R13_14C) - alam14 * R12_14C # deep ocean layer
    dR13_14C = ak2 * (R12_14C - (2* R13_14C) + R14_14C) - alam14 * R13_14C # deep ocean layer
    dR14_14C = ak2 * (R13_14C - (2* R14_14C) + R15_14C) - alam14 * R14_14C # deep ocean layer
    dR15_14C = ak2 * (R14_14C - (2* R15_14C) + R16_14C) - alam14 * R15_14C # deep ocean layer
    dR16_14C = ak2 * 2* (R15_14C - R16_14C) - alam14 * R16_14C  # deep ocean layer
    
    
    return(list(c(dRa_14C, dRb_14C,dRm_14C, dR1_14C, dR2_14C, dR3_14C, dR4_14C, dR5_14C,dR6_14C, dR7_14C, dR8_14C, dR9_14C, dR10_14C,
                  dR11_14C, dR12_14C, dR13_14C, dR14_14C, dR15_14C, dR16_14C), Q = Q, dRa_14C = dRa_14C, dRm_14C = dRm_14C
    )) 
  })
}

## Integration with 'ode': Use the "rk4" method if there are instabilities. Default intergrator is faster. 
# hmax is maximum time step - needs to be smaller than h1^2/(2 K) = ~ 1.8 (Courant criterion)
# 
# out <- ode(y = state, times = time, func = diffusion, parms = parameters,method = "rk4", hini=1.)
out <- ode(y = state, times = time, func = diffusion, parms = parameters,hmax=1.)

# convert 14Cratios into Delta-C14 (in permil; actually permil deviation from preindustrial steady state)
d14a = (out[,"Ra_14C"]/rs14) * 1000.
d14m = (out[,"Rm_14C"]/rs14) * 1000.
d14b = (out[,"Rb_14C"]/rs14) * 1000.
d14_1 = (out[,"R1_14C"]/rs14) * 1000.
d14_16 = (out[,"R16_14C"]/rs14) * 1000.

####### Plotting ######
par (mfrow=c(3, 2))
# 14C output
plot(out[,"time"], d14a,type="l", xlab="Year",ylab="Delta_14C_a",col=6)
plot(out[,"time"], d14b,type="l", xlab="Year",ylab="Delta_14C_b", col=6)
plot(out[,"time"], d14m,type="l", xlab="Year",ylab="Delta_14C_m", col=6)
plot(out[,"time"], d14_1,type="l", xlab="Year",ylab="Delta_14C_1", col=6)
plot(out[,"time"], d14_16,type="l", xlab="Year",ylab="Delta_14C_16", col=6)
plot(out[,"time"], out[,"Q"], type="l", col=2, ylab="Forcing")



