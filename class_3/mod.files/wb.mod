TITLE Wang-Buzsaki model

COMMENT

Sungho Hong, CNS Unit, Okinawa Inst Sci Tech. Mar 2008
ENDCOMMENT

UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
        (S) = (siemens)
}

? interface
NEURON {
        SUFFIX wb
        USEION na READ ena WRITE ina
        USEION k READ ek WRITE ik
        NONSPECIFIC_CURRENT il
        RANGE gnabar, gkbar, gl, el, gna, gk
        GLOBAL minf, hinf, ninf, mtau, htau, ntau
        THREADSAFE
}

PARAMETER {
        gnabar = .035 (S/cm2)	<0,1e9>
        gkbar = .009 (S/cm2)	<0,1e9>
        gl = .0001 (S/cm2)	<0,1e9>
        el = -65 (mV)
}

STATE {
        m h n
}

ASSIGNED {
        v (mV)
        celsius (degC)
        ena (mV)
        ek (mV)

        gna (S/cm2)
        gk (S/cm2)
        ina (mA/cm2)
        ik (mA/cm2)
        il (mA/cm2)
        minf hinf ninf
        mtau (ms) htau (ms) ntau (ms)
}

LOCAL mexp, hexp, nexp

? currents
BREAKPOINT {
        SOLVE states METHOD cnexp
        gna = gnabar*m*m*m*h
        ina = gna*(v - ena)
        gk = gkbar*n*n*n*n
        ik = gk*(v - ek)
        il = gl*(v - el)
}


INITIAL {
	rates(v)
	m = minf
	h = hinf
	n = ninf
}

? states
DERIVATIVE states {
        rates(v)
        m = minf
        h' = (hinf-h)/htau
        n' = (ninf-n)/ntau
}

LOCAL q10


? rates
PROCEDURE rates(v(mV)) { :Computes rate and other constants at current v.
                        :Call once from HOC to initialize inf at resting v.
        LOCAL  alpha, beta, sum
        TABLE minf, mtau, hinf, htau, ninf, ntau DEPEND celsius FROM -100 TO 100 WITH 200

UNITSOFF
        q10 = 5  :no temperature dependence
        :"m" sodium activation system
        alpha = .1 * vtrap(-(v+35),10)
        beta =  4 * exp(-(v+60)/18)
        sum = alpha + beta
        : mtau = 1/(q10*sum)
        minf = alpha/sum

        :"h" sodium inactivation system
        alpha = .07 * exp(-(v+58)/20)
        beta = 1 / (exp(-(v+28)/10) + 1)
        sum = alpha + beta
        htau = 1/(q10*sum)
        hinf = alpha/sum

        :"n" potassium activation system
        alpha = .01*vtrap(-(v+34),10)
        beta = .125*exp(-(v+44)/80)
        sum = alpha + beta
        ntau = 1/(q10*sum)
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
