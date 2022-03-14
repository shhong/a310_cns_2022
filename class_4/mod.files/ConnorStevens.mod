TITLE Connor-Stevens model

COMMENT
Connor-Stevens model adapted from Dayan and Abbott.
Implemented by Mari Shishikura, Kyoto University
March 2019
ENDCOMMENT

UNITS{
      (mA) = (milliamp)
      (mV) = (millivolt)
      (S) = (siemens)
}

NEURON{
      SUFFIX conste
      USEION na READ ena WRITE ina
      USEION k READ ek WRITE ik
      NONSPECIFIC_CURRENT il
      RANGE gnabar, gkbar, gl, el, gna, gk, ga, gabar
      GLOBAL minf, hinf, ninf, ainf, binf, mtau, htau, ntau, atau, btau
      THREADSAFE: assigned GLOBALs will be per Thread
}

PARAMETER {
      ena = 55 (mV)
      ek = -72 (mV)
      eashift = 3 (mV)
      gl = 3e-4 (S/cm2)
      el = -70 (mV)
      gnabar = .12 (S/cm2)
      gkbar = .02 (S/cm2)
      gabar = .0477 (S/cm2)
}

ASSIGNED {
      v (mV)
      gna (S/cm2)
      gk (S/cm2)
      ga (S/cm2)
      ina (mA/cm2)
      ik (mA/cm2)
      il (mA/cm2)
      minf hinf ninf ainf binf
      mtau (ms)
      htau (ms)
      ntau (ms)
      atau (ms)
      btau (ms)
}

STATE {
      m h n a b
}

BREAKPOINT{
      SOLVE states METHOD cnexp
      gna = gnabar*m*m*m*h
      ina = gna*(v-ena)

      gk = gkbar*n*n*n*n
      ga = gabar*a*a*a*b

      ik = gk*(v-ek) + ga*(v-ek+eashift)

      il = gl*(v-el)
}

INITIAL{
      rates(v)
      m = minf
      h = hinf
      n = ninf
      a = ainf
      b = binf
}

DERIVATIVE states {
      rates(v)
      m' = (minf-m)/mtau
      h' = (hinf-h)/htau
      n' = (ninf-n)/ntau
      a' = (ainf-a)/atau
      b' = (binf-b)/btau
}

PROCEDURE rates(v(mv)){
      LOCAL alpha, beta, sum, x
      TABLE minf, mtau, hinf, htau, ninf, ntau, ainf, atau, binf, btau FROM -100 TO 100 WITH 200
UNITSOFF
          :"m" sodium activation system
      alpha = .38*vtrap(-v-29.7, 10) :*(v+29.7)/(1-exp(-0.1*(v+29.7)))
      beta = 15.2*exp(-0.0556*(v+54.7))
      sum = alpha + beta
      mtau = 1/sum
      minf = alpha/sum
          :"h" sodium inactivation system
      alpha = 0.266*exp(-0.05*(v+48))
      beta = 3.8/(1+exp(-0.1*(v+18)))
      sum = alpha + beta
      htau = 1/sum
      hinf = alpha/sum
          :"n" potassium activation system
      alpha = .02*vtrap(-v-45.7, 10) :*(v+45.7)/(1-exp(-.1*(v+45.7)))
      beta = .25*exp(-.0125*(v+55.7))
      sum = alpha + beta
      ntau = 1/sum
      ninf = alpha/sum
          :"a" A-type potassium activation system
      x = (.0761*exp(.0314*(v+94.22)))/(1+exp(0.0346*(v+1.17)))
      ainf = pow(x, 1/3.0)
      atau = .3632 + 1.158/(1+exp(.0497*(v+55.96)))
          :"b" A-type potassium activaion system
      x = 1/(1+exp(.0688*(v+53.3)))
      binf = x*x*x*x
      btau = 1.24+2.678/(1+exp(.0624*(v+50)))
}


FUNCTION vtrap(x,y) {  :Traps for 0 in denominator of rate eqns.
        if (fabs(x/y) < 1e-6) {
                vtrap = y*(1 - x/y/2)
        }else{
                vtrap = x/(exp(x/y) - 1)
        }
}

UNITSON
