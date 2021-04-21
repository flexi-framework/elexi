import numpy as np
import matplotlib.pyplot as plt
from tikzplotlib import save as tikz_save

# variables
tauSGS   = 2*np.pi/32
kSGS     = 3/2
sigmaSGS = np.sqrt(2/3*kSGS)
betaSGS  = 1.
udiff    = [1,1,1]
urel     = udiff/np.sqrt(udiff[0]**2+udiff[1]**2+udiff[2]**2)

# tau parallel
tau = [tauSGS/np.sqrt(1.+betaSGS**2*(udiff[0]**2+udiff[1]**2+udiff[2]**2)/kSGS*3/2)]
# tau perpendicular
tau.append(tauSGS/np.sqrt(1.+4*betaSGS**2*(udiff[0]**2+udiff[1]**2+udiff[2]**2)/kSGS*3/2))

tend = 5.
t    = np.linspace(0.,tend,1000)
uSGS = np.zeros((3,len(t)))
w = np.zeros((3,3))
e = np.zeros((3,3))

for k in range(1,len(t)):
  dt = t[k]-t[k-1]
  for i in range(3):
    for j in range(3):
      if i==j:
        e[i,j] = np.exp(-dt/tau[1])+(np.exp(-dt/tau[0])-np.exp(-dt/tau[1]))*urel[i]*urel[j]
        w[i,j] = sigmaSGS*np.sqrt(1-np.exp(-2*dt/tau[1]))+(sigmaSGS*np.sqrt(1-np.exp(-2*dt/tau[0]))-sigmaSGS*np.sqrt(1-np.exp(-2*dt/tau[1])))*urel[i]*urel[j]
      else:
        e[i,j] = (np.exp(-dt/tau[0])-np.exp(-dt/tau[1]))*urel[i]*urel[j]
        w[i,j] = (sigmaSGS*np.sqrt(1-np.exp(-2*dt/tau[0]))-sigmaSGS*np.sqrt(1-np.exp(-2*dt/tau[1])))*urel[i]*urel[j]

for k in range(1,len(t)):
  for j in range(3):
    uSGS[:,k] = uSGS[:,k]+e[:,j]*uSGS[j,k-1]+w[:,j]*np.random.normal(0,1)

plt.figure(figsize=(8,6),dpi=300)
plt.xlim(0,5)
plt.ylim(0,8)
plt.xlabel('time [s]',usetex=True)
plt.ylabel("u'(t)"   ,usetex=True)

plt.plot(t,uSGS[0,:])
plt.show()
#tikz_save('SGS_Analytic_Velo.tikz')
plt.close()

tmp = 0
dt  = tend/1000
steps = int((tend-0)/dt)
y = np.zeros((steps))
for k in range(1,len(t)):
    dt   = t[k]-t[k-1]
    tmp  = y[k-1] + uSGS[0,k]*dt
    y[k] = tmp

plt.figure(figsize=(8,6),dpi=300)
plt.xlim(0,5)
plt.ylim(0,5)
plt.xlabel('time [s]'   ,usetex=True)
plt.ylabel("Path length",usetex=True)
plt.plot(t,y)
plt.show()
#tikz_save('SGS_Analytic_Path.tikz')
plt.close()
