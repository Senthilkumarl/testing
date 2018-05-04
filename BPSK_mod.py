
# coding: utf-8
__author__ = 'Senthilkumar Lakshmanan'
__email__ = 'senthilkumarl@live.com'
__date__ = '11-04-2018'

########################BPSK modulation########################################


from math import sqrt, pi
import numpy as np
import scipy
import matplotlib.pyplot as plt


def add_AWGN(x,nPr):
    return x + sqrt(nPr)*np.random.randn(len(x))

def dot_prod(x, Phi):
    return sum([x[i] * Phi[i] for i in range(len(x))])
def vec_prod(x, Phi):
    return [x[i] * Phi[i] for i in range(len(x))]

def bpskMod(params,bit,t):
    f = params["f"]
    #f=10
    T = 1/f
    Phi1 = params["Phi"][0](t)#A*np.cos( 2.0*pi*f*t + pi/4)
    Phi2 = params["Phi"][1](t)#A*np.cos( 2.0*pi*f*t - pi/4)
    # B = 2/sqrt(T)
    # x0 = B*np.sin( 2.0*pi*f*t)
    # x1 = -B*np.sin( 2.0*pi*f*t)
    # x0 = Phi1 - Phi2
    # x1 = Phi2 - Phi1

    Phi = np.array([Phi1, Phi2])
    if bit == 1:
        x = params["sym"][0]#np.array([1, -1])
    else:
        x = params["sym"][1]#np.array([-1, 1])
    x_t = dot_prod(x, Phi)
    # x = np.array([-1, 1])
    # x1 = dot_prod(x, Phi)
    nPr = params["nPr"]
    x_t = add_AWGN(x_t,nPr)
    if params["plt"]:
        plt.subplot(211)
        plt.plot(t,x_t)
    return x_t
    #plt.plot(t,Phi1,t,Phi2)

def detector(detI,detQ,params):
    sym = params["sym"]
    symI = [sym[i][0] for i in range(len(sym))]
    symQ = [sym[i][1] for i in range(len(sym))]
    #print(symI)
    UBI = max(symI)
    LBI = min(symI)
    UBQ = max(symQ)
    LBQ = min(symQ)
    if detI > 0:
        I = UBI #1
    else:
        I = LBI #-1

    if detQ > 0:                        #"sym":[np.array([1, -1]),np.array([-1, 1])]
        Q = UBQ #1
    else:
        Q = LBQ #-1

    if I == params["sym"][1][0] and Q == params["sym"][1][1]:#I == -1 and Q == 1:
        return 0
    else:
        return 1

def bpsdDemod(x,t,params):
    f = params["f"]
    dt = params["dt"]
    T = 1/f
    A=sqrt(2/T)
    Phi1 = params["Phi"][0](t)#A*np.cos( 2.0*pi*f*t + pi/4)
    Phi2 = params["Phi"][1](t)#A*np.cos( 2.0*pi*f*t - pi/4)
    yI = vec_prod(np.asarray(x),Phi1)
    detI = np.trapz(yI,dx=dt)#/abs(np.trapz(Phi1*Phi1))
    yQ = vec_prod(np.asarray(x),Phi2)
    detQ = np.trapz(yQ,dx=dt)#/np.trapz(Phi2*Phi2)
    plt.subplot(212)
    plt.plot(detI,detQ,"*g")
    return detector(detI,detQ,params)

def plot_canstn(sym):
    symI = [sym[i][0] for i in range(len(sym))]
    symQ = [sym[i][1] for i in range(len(sym))]
    #print(sym,symI,symQ)
    plt.subplot(212)
    plt.plot(symI,symQ,"*r")

def main():
    dt = 0.001
    f = 10
    T = 1/f
    Tb = 1*T
    bit = 0
    N = 10
    Pr = 1

    A=sqrt(2/T)*Pr
    params = {
            "f":f,
            "Tb":Tb,
            "dt":dt,
            "N":N,
            "Phi":[lambda t:A*np.cos( 2.0*pi*f*t + pi/4), lambda t:A*np.cos( 2.0*pi*f*t - pi/4)],#[lambda t:A*np.cos( 2.0*pi*f*t), lambda t:A*np.cos( 2.0*pi*f*t+pi/2)],#
            "sym":[np.array([1, -1]),np.array([-1, 1])],#[np.array([-1, 0]),np.array([1, 0])],#[np.array([1, 1]),np.array([-1, -1])],#
            "nPr":1,
            "ax":[],
            "plt":True,

            }

    bits =np.random.randint(2, size=N)# [0,0,0,1,1,0,1,1,0,1]
    det_bits = []
    for i in range(N):
        t=np.arange(i*Tb,(i+1)*Tb+dt,dt)
        x_t = bpskMod(params,bits[i],t)
        det_bit = bpsdDemod(x_t,t,params)
        det_bits.append(det_bit)

    BER = sum( [0 if bits[i] == det_bits[i] else 1 for i in range(N)] ) / len(bits)
    print(BER)

    print(det_bits)


    if params["plt"]:
        plt.subplot(211)
        plt.title("$T_b=$"+str(Tb)+", T="+str(T))
        plot_canstn(params["sym"])
        plt.show()

if __name__ == '__main__':
    main()
