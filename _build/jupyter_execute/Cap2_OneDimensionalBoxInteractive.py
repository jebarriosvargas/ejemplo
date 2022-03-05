#!/usr/bin/env python
# coding: utf-8

# # Partícula en una caja -- analítica

# In[1]:


from pylab import *
import seaborn as sns
from ipywidgets import *


# \begin{align}
# E_n = \bigg( \frac{h^2}{8m}\bigg) \frac{1}{\ell^2} n^2
# \end{align}
# 
# \begin{align}
# \psi(x) = \sqrt{\frac{2}{\ell}} \sin\bigg( \frac{n\pi x}{\ell} \bigg)
# \end{align}

# In[2]:


π          = pi
# Physical quantities
e          = 1.602e-19
h          = 4.135667               # eV fs
c          = 299792458              # [m/s]
massfactor = e/c/c                  # 1 eV/c^2 = 1.79e-36 kg
me         = 9.109e-31/massfactor   # [eV/c^2] = 0.5x10^6 eV/c^2   
c_nmfs     = 299.792458             # [nm/fs]
Eλ         = h*h*c_nmfs*c_nmfs/8/me # eV nm^2

# Number of levels to show
nmax       = 5
n          = arange(nmax)+1

ℓ          = 1                      # nm

E          = Eλ/ℓ/ℓ*n*n

x          = linspace(0,ℓ,1000)
ψ          = array( [ sqrt(2/ℓ)*sin(nv*π*x/ℓ) for nv in n] )
ψ2         = ψ*ψ


# In[3]:


fig,axes = plt.subplots(ncols=2,nrows=1,figsize=(11,5))

axes[0].set_xlim(-.4*ℓ  ,1.4*ℓ)

color_list = sns.cubehelix_palette(nmax)

for i in arange(nmax):
    axes[0].plot(x,E[i]+ψ[i],label=r"$E_{0}$={1:>8.3f}".format(i+1,E[i]),c=color_list[i])

for i in arange(nmax):
    axes[0].plot([0 ,ℓ],[E[i],E[i]],'--',c=color_list[i])

    

axes[0].plot([0,0],[0,(ψ.max()+E[-1])*1.1],c="Gray")
axes[0].plot([ℓ,ℓ],[0,(ψ.max()+E[-1])*1.1],c="Gray")
axes[0].plot([0,ℓ],[0,0],c="Gray")

axes[0].set_xlabel("$x$ [nm]",fontsize=16)
axes[0].set_ylabel("$E$ [eV]",fontsize=16)
axes[0].legend(loc=1)
axes[0].set_title(r"$\psi$")

axes[1].set_xlim(-.4*ℓ  ,1.4*ℓ)

color_list = sns.cubehelix_palette(nmax)

for i in arange(nmax):
    axes[1].plot(x,E[i]+ψ2[i],label=r"$E_{0}$={1:>8.3f}".format(i+1,E[i]),c=color_list[i])

for i in arange(nmax):
    axes[1].plot([0 ,ℓ],[E[i],E[i]],'--',c=color_list[i])

    

axes[1].plot([0,0],[0,(ψ2.max()+E[-1])*1.1],c="Gray")
axes[1].plot([ℓ,ℓ],[0,(ψ2.max()+E[-1])*1.1],c="Gray")
axes[1].plot([0,ℓ],[0,0],c="Gray")

axes[1].set_xlabel("$x$ [nm]",fontsize=16)
axes[1].set_ylabel("$E$ [eV]",fontsize=16)
axes[1].legend(loc=1)
axes[1].set_title(r"$|\psi|^2$")


fig.tight_layout()
# fig.savefig("Cap2_OneDimensionalBox.pdf")


# ## Evolución temporal

# In[4]:


π          = pi
# Physical quantities
e          = 1.602e-19
h          = 4.135667               # eV fs
ħ          = h/2/π                  # eV fs
c          = 299792458              # [m/s]
massfactor = e/c/c                  # 1 eV/c^2 = 1.79e-36 kg
me         = 9.109e-31/massfactor   # [eV/c^2] = 0.5x10^6 eV/c^2   
c_nmfs     = 299.792458             # [nm/fs]
Eλ         = h*h*c_nmfs*c_nmfs/8/me # eV nm^2

# Number of levels to show
nmax       = 5
n          = arange(nmax)+1

ℓ          = 1                      # nm

E          = Eλ/ℓ/ℓ*n*n

x          = linspace(0,ℓ,1000)
ψ          = array( [ sqrt(2/ℓ)*sin(nv*π*x/ℓ) for nv in n] )


# In[5]:


color_list = sns.cubehelix_palette(nmax)


# In[6]:


@interact(t=(0,10,0.01))
def Evolution(t=0.0):
    fig,axes = plt.subplots(ncols=3,nrows=1,figsize=(16,5))
    axes[0].set_xlim(-.4*ℓ  ,1.4*ℓ)
    axes[1].set_xlim(-.4*ℓ  ,1.4*ℓ)
    axes[2].set_xlim(-.4*ℓ  ,1.4*ℓ)
    
    axes[0].set_ylim(-2,2*nmax+2)
    axes[1].set_ylim(-2,2*nmax+2)
    axes[2].set_ylim(-2,2*nmax+2)
    

    
    for nv in n:
        
        i = nv - 1
        E  = Eλ/ℓ/ℓ*nv*nv
        Ψe = exp(-1.J*E*t/ħ)*sqrt(2/ℓ)*sin(nv*π*x/ℓ)
        Ψe2 = conjugate(Ψe)*Ψe
        
        axes[0].plot(x,E+Ψe.real,label=r"$E_{0}$={1:>8.3f}".format(i+1,E),c=color_list[i])
        axes[0].plot([0 ,ℓ],[E,E],'--',c=color_list[i])
        
        axes[1].plot(x,E+Ψe.imag,label=r"$E_{0}$={1:>8.3f}".format(i+1,E),c=color_list[i])
        axes[1].plot([0 ,ℓ],[E,E],'--',c=color_list[i])
        
        axes[2].plot(x,E+Ψe2.real,label=r"$E_{0}$={1:>8.3f}".format(i+1,E),c=color_list[i])
        axes[2].plot([0 ,ℓ],[E,E],'--',c=color_list[i])
        
    axes[0].plot([0,0],[0,2*nmax+2],c="Gray")
    axes[0].plot([ℓ,ℓ],[0,2*nmax+2],c="Gray")
    axes[0].plot([0,ℓ],[0,0],c="Gray")
    
    axes[0].set_xlabel("$x$ [nm]",fontsize=16)
    axes[0].set_ylabel("$E$ [eV]",fontsize=16)
    axes[0].legend(loc=1)
    axes[0].set_title(r"Re $(\psi)$")
    
    axes[1].plot([0,0],[0,2*nmax+2],c="Gray")
    axes[1].plot([ℓ,ℓ],[0,2*nmax+2],c="Gray")
    axes[1].plot([0,ℓ],[0,0],c="Gray")
    
    axes[1].set_xlabel("$x$ [nm]",fontsize=16)
    axes[1].set_ylabel("$E$ [eV]",fontsize=16)
    axes[1].legend(loc=1)
    axes[1].set_title(r"Im $(\psi)$")
    
    axes[2].plot([0,0],[0,2*nmax+2],c="Gray")
    axes[2].plot([ℓ,ℓ],[0,2*nmax+2],c="Gray")
    axes[2].plot([0,ℓ],[0,0],c="Gray")
    
    axes[2].set_xlabel("$x$ [nm]",fontsize=16)
    axes[2].set_ylabel("$E$ [eV]",fontsize=16)
    axes[2].legend(loc=1)
    axes[2].set_title(r"$|\psi|^2$")
    return
        


# ## Ortogonalidad

# In[7]:


ψ  = array( [ sqrt(2/ℓ)*sin(nv*π*x/ℓ) for nv in n] )

print( np.trapz( np.conj(ψ[0])*ψ[0],x ) )


# ## Operadores

# \begin{align*}
# \langle x \rangle &=&\frac{a}{2} \\
# \langle x^2 \rangle &=& \bigg( \frac{1}{3}-\frac{1}{2\pi^2} \bigg)a^2
# \end{align*}

# In[8]:


xmean = np.trapz( np.conj(ψ[0])*x*x*ψ[0],x  )
print (xmean)
print (1/3-1/2/pi/pi)


# \begin{align*}
# \langle p_x \rangle &=&0 \\
# \langle p_x^2 \rangle &=& \bigg( \frac{h^2}{4a^2} \bigg)
# \end{align*}

# In[9]:


ħ = h/2/π
dψ = np.gradient(ψ[0],x,edge_order=2)
pmean = np.trapz( np.conj(ψ[0])*dψ,x )
print (pmean)


# In[10]:


d2ψ = np.gradient(dψ,x,edge_order=2)
p2mean = np.trapz( np.conj(ψ[0])*d2ψ,x )
print (-ħ*ħ*p2mean)


# In[11]:


h*h/4


# In[12]:


plot(x,ψ[0])
plot(x,dψ)
plot (x,d2ψ)


# In[ ]:





# In[ ]:




