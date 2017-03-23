from numpy import exp
from matplotlib import pyplot as plt

class Ponto(object):
    def __init__(self,L,OD,Q,T):
        self.L = L
        self.OD = OD
        self.Q = Q
        self.T = T
    def __add__(self,other):
        Q = self.Q+other.Q
        T = (self.Q*self.T+other.Q*other.T)/Q
        OD = (self.Q*self.OD+other.Q*other.OD)/Q
        L = (self.Q*self.L+other.Q*other.L)/Q
        return Ponto(L,OD,Q,T)
class Trecho(object):
    def __init__(self,x,montante,h,A,Kd=0.35):
        self.x = x
        self.montante = montante
        self.h = h
        self.A = A
        self.U = (montante.Q/A)
        self.t = (1000/self.U)/86400
        self.Kd = Kd
        self.Kdt = Kd*1.047**(montante.T-20)
        self.Ka = 3.93*((self.U)**(1/2)/((h)**(1.5)))
        self.Kat = self.Ka*(1.024)**(self.montante.T-20)
        self.Ks = 0
        self.Kr = self.Kdt + self.Ks
        self.ODsat = exp(-139.34+(157570.1/(montante.T+273))-(66423080/((montante.T+273)**2))+(12438000000/((montante.T+273)**3))-(862196400000/((montante.T+273)**4)))
    def por_km(self):
        ODs = [self.montante.OD,]
        Ls = [self.montante.L,]
        Ds = [self.ODsat - self.montante.OD,]
        for i in range(self.x):
            OD_ant = ODs[i]
            L_ant = Ls[i]
            D = Ds[i]*exp(-self.Ka*self.t)+(self.Kdt*L_ant/(self.Ka-self.Kr))*(exp(-self.Kr*self.t)-exp(-self.Ka*self.t))
            Ds.append(D)
            ODs.append(self.ODsat-D)
            Ls.append(L_ant*exp((-self.Kr)*self.t))
        self.jusante = Ponto(Ls[-1],ODs[-1],self.montante.Q,self.montante.T)
        return self.jusante,Ls,Ds,ODs  
def unir_listas(L,A1,A2):
    L.extend(A1)
    L.extend(A2)
#Criando pontos de contribuição e trechos
P1 = Ponto(2,7.5,4.21,20.1)
C1 = Ponto(350,0,.35,28)
M1 = P1+C1
T1 = Trecho(20,M1,1.15,14.10)
J1,Ls,Ds,ODs = T1.por_km()
C2 = Ponto(500,0,.5,30)
M2 = C2 + J1
T2 = Trecho(20,M2,1.15,14.10)
J2,Ls2,Ds2,ODs2 = T2.por_km()
C3 = Ponto(5,0,.5,30)
M3 = J2 + C3
T3 = Trecho(60,M3,1.35,17.70)
J3,Ls3,Ds3,ODs3 = T3.por_km()
unir_listas(Ls,Ls2,Ls3)
unir_listas(Ds,Ds2,Ds3)
unir_listas(ODs,ODs2,ODs3)
#cria o gráfico
x = range(len(Ls),0,-1)
plt.plot(x,Ls,c = "k",label = "DBO")
plt.plot(x,ODs,c = "b",label = "OD")
#plt.plot(x,Ds,c = "r",label = "Déficit")
plt.legend(loc="best")
ax = plt.gca()
ax.invert_xaxis()
plt.xlabel("Distância (m)")
plt.ylabel("Concentração (mg/L)")
plt.show()

