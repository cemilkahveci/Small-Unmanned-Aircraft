import math 
import numpy as np

class Aerosonde():
    def __init__(self):
        self.m = 13.5 #kg

        Jx = 0.844 #kgm^2
        Jy = 1.135
        Jz = 1.759
        Jxz = 0.1204

        self.Jx = Jx
        self.Jy = Jy
        self.Jz = Jz
        self.Jxz = Jxz

        self.S = 0.55 #m^2
        self.S_prop = 0.2027
        self.b = 2.8956 #m
        self.c = 0.18994 
        self.rho = 1.2682 #kg/m^3
        self.k_motor = 80
        self.k_Tp = 0
        self.k_ohm = 0
        self.e = 0.9
        self.AR = self.b**2/self.S
        self.g = 9.81

        #longitudinal coefficients
        self.c_L0 = 0.28
        self.c_D0 = 0.03
        self.c_M0 = -0.02338
        self.c_La = 3.45
        self.c_Da = 0.3
        self.c_Ma = -0.38
        self.c_Lq = 0
        self.c_Dq = 0
        self.c_Mq = -3.6
        self.c_Lde = -0.36
        self.c_Dde = 0
        self.c_Mde = -0.5
        self.c_prop = 1
        self.M = 50
        self.a0 = 0.4712
        self.epsilon = 0.1592
        self.c_Dp = 0.0437
        self.c_ndr = -0.032

        #lateral coefficients
        self.c_Y0 = 0
        self.c_l0 = 0
        self.c_n0 = 0
        self.c_Yb = -0.98
        self.c_lb = -0.12
        self.c_nb = 0.25
        self.c_Yp = 0
        self.c_lp = -0.25
        self.c_np = 0.022
        self.c_Yr = 0
        self.c_lr = 0.14
        self.c_nr = -0.35
        self.c_Yda = 0
        self.c_lda = 0.08
        self.c_nda = 0.06
        self.c_Ydr = -0.17
        self.c_ldr = 0.105

        #Gamma for 6dof
        G = Jx*Jz - Jxz**2
        self.G1 = (Jxz*(Jx-Jy+Jz))/G
        self.G2 = (Jz*(Jz-Jy) + Jxz**2)/G
        self.G3 = Jz/G
        self.G4 = Jxz/G
        self.G5 = (Jz-Jx)/Jy
        self.G6 = Jxz/Jy
        self.G7 = (Jx*(Jx-Jy) + Jxz**2)/G
        self.G8 = Jx/G

        self.c_p0 = self.G3*self.c_l0 + self.G4*self.c_n0
        self.c_pbeta = self.G3*self.c_lb + self.G4*self.c_nb
        self.c_pp = self.G3*self.c_lp + self.G4*self.c_np
        self.c_pr = self.G3*self.c_lr + self.G4*self.c_nr
        self.c_pda = self.G3*self.c_lda + self.G4*self.c_nda
        self.c_pdr = self.G3*self.c_ldr + self.G4*self.c_ndr
        self.c_r0 = self.G4*self.c_l0 + self.G8*self.c_n0
        self.c_rbeta = self.G4*self.c_lb + self.G8*self.c_nb
        self.c_rp = self.G4*self.c_lp + self.G8*self.c_np
        self.c_rr = self.G4*self.c_lr + self.G8*self.c_nr
        self.c_rda = self.G4*self.c_lda + self.G8*self.c_nda
        self.c_rdr = self.G4*self.c_ldr + self.G8*self.c_ndr    

    def Derivatives(self):
        self.Va = math.sqrt(self.u**2 + self.v**2 + self.w**2)
        self.alpha = math.atan(self.w/self.u)
        self.beta = math.asin(self.v/self.Va)

        sigma_a = (1 + np.exp(-self.M*(self.alpha - self.a0)) + np.exp(self.M*(self.alpha + self.a0)))/((1 + np.exp(-self.M*(self.alpha - self.a0)))*np.exp(self.M*(self.alpha + self.a0)))
        func_c_La = (1-sigma_a)*(self.c_L0+self.c_La*self.alpha) + sigma_a * (2*np.sign(self.alpha)*np.sin(self.alpha)**2*np.cos(self.alpha))
        func_c_Da = self.c_Dp + (self.c_L0 + self.c_La*self.alpha)**2/(np.pi*self.e*self.AR)
        c_Xa = -func_c_Da*np.cos(self.alpha) + func_c_La*np.sin(self.alpha)
        c_Xq = -self.c_Dq*np.cos(self.alpha) + self.c_Lq*np.sin(self.alpha)
        c_Xde = -self.c_Dde*np.cos(self.alpha) + self.c_Lde*np.sin(self.alpha)
        c_Z = -func_c_Da*np.sin(self.alpha) - func_c_La*np.sin(self.alpha)
        c_Zq = -self.c_Dq*np.sin(self.alpha) - self.c_Lq*np.sin(self.alpha)
        c_Zde = -self.c_Dde*np.sin(self.alpha) - self.c_Lde*np.sin(self.alpha)        



        self.dot_pn = (self.e1**2+self.e0**2-self.e2**2-self.e3**2)*self.u + 2*(self.e1*self.e2-self.e3*self.e0)*self.v + 2*(self.e1*self.e3+self.e2*self.e0)*self.w
        self.dot_pe = 2*(self.e1*self.e2+self.e3*self.e0)*self.u + (self.e2**2+self.e0**2-self.e1**2-self.e3**2)*self.v + 2*(self.e2*self.e3-self.e1*self.e0)*self.w
        self.dot_pd = 2*(self.e1*self.e3-self.e2*self.e0)*self.u + 2*(self.e2*self.e3+self.e1*self.e0)*self.v + (self.e3**2+self.e0**2-self.e2**2-self.e1**2)*self.w

        self.dot_u = self.r*self.v - self.q*self.w + 2*self.g*(self.e1*self.e3-self.e2*self.e0) + ((self.rho*(self.Va**2)*self.S)/(2*self.m))*(c_Xa + c_Xq*(self.c*self.q)/(2*self.Va) + c_Xde*self.delta_e) + ((self.rho*self.S_prop*self.c_prop)/(2*self.m))*((self.k_motor*self.delta_t)**2 - self.Va**2)
        self.dot_v = self.p*self.w - self.r*self.u + 2*self.g*(self.e2*self.e3 + self.e1*self.e0) + ((self.rho*(self.Va**2)*self.S)/(2*self.m))*(self.c_Y0 + self.c_Yb*self.beta + (self.c_Yp*self.b*self.p)/(2*self.Va) + (self.c_Yr*self.b*self.r)/(2*self.Va) + self.c_Yda*self.delta_a + self.c_Ydr*self.delta_r)
        self.dot_w = self.q*self.u - self.p*self.v + self.g(self.e3**2 + self.e0**2 - self.e1**2 -self.e2**2) + ((self.rho*(self.Va**2)*self.S)/(2*self.m))*(c_Z + (c_Zq*self.c*self.q)/(2*self.Va) + c_Zde*self.delta_e)

        self.dot_e0 = .5*(-self.p*self.e1 -self.q*self.e2 - self.r*self.e3)
        self.dot_e1 = .5*(self.p*self.e0 + self.r*self.e2 - self.q*self.e3)
        self.dot_e2 = .5*(self.q*self.e0 - self.r*self.e1 + self.p*self.e3)
        self.dot_e3 = .5*(self.r*self.e0 + self.q*self.e1 - self.p*self.e2)

        self.dot_p = self.G1*self.p*self.q - self.G2*self.q*self.r + (0.5*self.rho*self.Va*self.Va*self.S*self.b)*(self.c_p0 + self.c_pbeta*self.beta + (self.c_pp*self.b*self.p)/(2*self.Va) + (self.c_pr*self.b*self.r)/(2*self.Va) + self.c_pda*self.delta_a + self.c_pdr*self.delta_r)
        self.dot_q = self.G5*self.p*self.r - self.G6*(self.p**2-self.r**2) + ((self.rho*self.Va*self.Va*self.S*self.c)/(2*self.Jy))*(self.c_M0 + self.c_Ma*self.alpha + (self.c_Mq*self.c*self.q)/(2*self.Va) + self.c_Mde*self.delta_e)
        self.dot_r = self.G7*self.p*self.q - self.G1*self.q*self.r + (0.5*self.rho*self.Va*self.Va*self.S*self.b)*(self.c_r0 + self.c_rbeta*self.beta + (self.c_rp*self.b*self.p)/(2*self.Va) + (self.c_rr*self.b*self.r)/(2*self.Va) + self.c_rda*self.delta_a + self.c_rdr*self.delta_r)
