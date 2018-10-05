
# coding: utf-8

# In[1]:

import numpy as np
from scipy import optimize

def f4(self,dx):
    return (fun(self+dx/2)-fun(self-dx/2))/dx
def f2(self,dx):
    return 2*(fun(self+dx/4)-fun(self-dx/4))/dx
def fun(self):
    return 1+self**2-self**3
def df(self):
    return 2*self-3*self**2

class Derivada:
    
    def __init__(self,f,metodo='adelante',dx=0.001):
        self.f=f
        self.metodo= metodo
        self.dx=dx
    def calc(self,x):
        self.x=x
        if self.metodo=='adelante':
            while self.dx>1e-7:
                f=(fun(x+self.dx)-fun(x))/self.dx
                self.dx=self.dx/2
            return f
        elif self.metodo=='extrapolada':
            while self.dx>1e-7:
                f=(4*f4(x,self.dx)-f2(x,self.dx))/3
                self.dx=self.dx/2
            return f
        
        elif self.metodo=='segunda':
            while self.dx>1e-7:
                f=(fun(x+self.dx)+fun(x-self.dx)-2*fun(x))/self.dx**2
                self.dx=self.dx/2
                return f
        elif self.metodo=='central':
            while self.dx>1e-7:
                f=(fun(x+0.5*self.dx)-fun(x-0.5*self.dx))/self.dx
                self.dx=self.dx/2
                return f
        else:
            return "Metodo no ha sido implementado"

class Zeros:
    def __init__(self,f,metodo,error=1e-4,max_iter=100):
        self.f=fun
        self.metodo=metodo
        self.error=error
        self.max_iter=max_iter
    def zero(self,vi):
        if self.metodo=='bisectriz':
            self.vi=vi
            vil=[]
            for k in range(len(vi)):
                vil.append(self.vi[k])
            c=vil[0]
            #rint(c)
            for i in range(self.max_iter): 
                c = (vil[0]+vil[1])/2
                
                if (fun(c)*fun(vil[0]) < 0): 
                    vil[1]=c 
                else:
                    vil[0] =c   
                if abs(fun(vil[0])<=self.error): 
                    break
            return c
        elif self.metodo=='newton':
            self.vi=vi
            vil=vi
            for i in range(self.max_iter):
                vil=vil-fun(vil)/df(vil)
                if abs(fun(vil))<self.error:
                    break
            return vil
        elif self.metodo=='interpolacion':
            self.vi=vi
            vil=[]
            for k in range(len(vi)):
                vil.append(self.vi[k])
            for i in range(self.max_iter):
                x=(vil[0]*fun(vil[1])-vil[1]*fun(vil[0]))/(fun(vil[1])-fun(vil[0]))
                if abs(fun(x))<self.error:
                    break
                elif fun(vil[0])*fun(x)<0:
                    vil[1]=x
                else:
                    vil[0]=x
            
            return x
        elif self.metodo=='newton-sp':
            self.vi=vi
            x0=vi
            args=()
            tol=self.error
            maxiter=self.max_iter
            r = optimize.newton(fun, x0, df,args,tol,maxiter)
            return r
        elif self.metodo=='brentq-sp':
            self.vi=vi
            vil=[]
            for i in range(len(vi)):
                vil.append(vi[i])
            a1=vil[0]
            b1=vil[1]
            args=()
            r=optimize.brentq(fun,a1,b1,args,xtol=self.error,maxiter=self.max_iter)
            return r
        elif self.metodo=='fsolve-sp':
            self.vi=vi
            x0=vi
            args=()
            r=optimize.fsolve(fun,x0,args,xtol=self.error,factor=self.max_iter)
            return str(r)
            
        
    def __getitem__(self, key):
        return self.Zeros[key]
            
            

                
        
if __name__ == "__main__":
    #Calculo de la derivada y raiz de la funcion -x^3+x^2+1 donde se tomo la derivada evaluada en 1 y la raiz cerca a 1
    f=0
    dev=Derivada(f,'adelante',dx=0.001).calc(1)
    print(dev)
    zep=Zeros(f,'brentq-sp',error=1e-4,max_iter=100).zero((1,2))
    print(zep)


# In[ ]:



