TITLE mease.mod

COMMENT
ENDCOMMENT

UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
        (S) = (siemens)
}

NEURON{
        SUFFIX mease
        USEION na READ ena WRITE ina
        USEION k READ ek WRITE ik
        NONSPECIFIC_CURRENT il
        RANGE gnabar, gkbar, gl, el, gna, gk
        GLOBAL minf, hinf, ninf, mtau, htau, ntau
        THREADSAFE : assigned GLOBALs will be per thread
}

PARAMETER {
      ena = 50 (mV)
      ek = -77 (mV)
      gl = .000025 (S/cm2)
      el = -70 (mV)
      gnabar = 14e-2 (S/cm2)
      gkbar = 14e-2 (S/cm2)
}

ASSIGNED {
      v (mV)

      gna (S/cm2)
      gk(S/cm2)
      ina (mA/cm2)
      ik(mA/cm2)
      il(mA/cm2)
      minf hinf ninf
      mtau (ms)
      htau (ms)
      ntau(ms)
}

STATE {
      m h n
}

BREAKPOINT {
      SOLVE states METHOD cnexp
      gna = gnabar*m*m*m*h
      ina = gna*(v - ena)
      gk = gkbar*n
      ik = gk*(v-ek)
      il = gl*(v-el)
}

INITIAL {
  rates(v)
  m = minf
  h = hinf
  n = ninf
}

DERIVATIVE states {
    rates(v)
    m' = (minf-m)/mtau
    h' = (hinf-h)/htau
    n' = (ninf-n)/ntau
}

PROCEDURE rates(v(mv)){  :Computes rate and other constants at current v.
                        :Call once from HOC to initialize inf at rensting v.
      LOCAL alpha, beta, sum
      TABLE minf, mtau, hinf, htau, ninf, ntau FROM -100 TO 100 WITH 200
UNITSOFF
          :"m" sodium activation system
      alpha = .182*vtrap(-v-35, 9) : (v+35)/(1-exp(-(v+35)/9))
      beta =  .124*vtrap(v+35, 9) : (v+35)/(1-exp((v+35)/9))
      sum = alpha + beta
      mtau = 1/sum
      minf = alpha/sum
          :"h" sodium inactivation system
      alpha = .024*vtrap(-v-50, 5) :(v+50)/(1-exp(-(v+50)/5))
      beta =  .0091*vtrap(v+75, 5) :(v+75)/(1-exp((v+75)/5))
      sum = alpha + beta
      htau = 1/sum
      hinf=1/(1+exp((v+65)/6.2))
         :"n" potassium activation system
      alpha = .020*vtrap(-v+20, 9) :(v-20)/(1-exp(-(v-20)/9))
      beta =  .002*vtrap(v-20, 9)  :(v-20)/(1-exp((v-20)/9)))
      sum = alpha + beta
      ntau = 1/sum
      ninf = alpha/sum
}

FUNCTION vtrap(x,y) {  :Traps for 0 in denominator of rate eqns.
        if (fabs(x/y) < 1e-6) {
                vtrap = y*(1 - x/y/2)
        }else{
                vtrap = x/(exp(x/y) - 1)
        }
}

UNITSON
