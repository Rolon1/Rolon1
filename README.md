import numpy as np
import matplotlib.pyplot as plt
import math as m
def f(x):
	return 2*m.sin(0.25*x)+m.cos(4*x)
a=2
b=4
n=3
h=(b-a)/n
xmatrix=np.zeros(n+1)
ymatrix=np.zeros(n+1)
for i in range(int(n)+1):
	xmatrix[i]=a+i*h
	ymatrix[i]=f(xmatrix[i])
print(xmatrix)
print(ymatrix)
def creatematrix(n):
	C=np.zeros((n+1,n+1))
	for i in range(n+1):
		for j in range(n+1):
			C[i,j]=xmatrix[i]**j
	return C
C=creatematrix(n)
print(C)
A=np.linalg.solve(C,ymatrix)
print(A)
print(C[1,1])
def Lagrandg(z):
	L=0
	a=np.zeros(n+1)
	b=np.zeros(n+1)
	for i in range(n+1):
		a[i]=1
		b[i]=1
		for j in range(n+1):
			if i!=j:
				a[i]=a[i]*(z - xmatrix[j])
				b[i]=b[i]*(xmatrix[i] - xmatrix[j])
		L=L+ymatrix[i]*a[i]/b[i]
	return L
def Newton(t):
	dy0=ymatrix[1]-ymatrix[0]
	d2y0=ymatrix[2]-2*ymatrix[1]+ymatrix[0]
	d3y0=ymatrix[3]-3*ymatrix[2]+3*ymatrix[1]-ymatrix[0]	
	return ymatrix[0]+(t-xmatrix[0])*dy0/h+(t-xmatrix[0])*(t-xmatrix[1])*d2y0/(2*h**2)+(t-xmatrix[0])*(t-xmatrix[1])*(t-xmatrix[2])*d3y0/(6*h**3)
def Eps(t):
	return abs(1-Lagrandg(t)/f(t))*100
x=np.linspace(2, 4, 50)
y1=np.linspace(2, 4, 50)
y2=[A[0]+A[1]*i+A[2]*i**2+A[3]*i**3 for i in x]
y3=np.linspace(2, 4, 50)
y4=np.linspace(2, 4, 50)
y5=np.linspace(2, 4, 50)
for i in range(50):
	y1[i]=Lagrandg(x[i])
	y3[i]=f(x[i])
	y4[i]=Newton(x[i])
	y5[i]=Eps(x[i])
plt.plot(x,y1, label='L')
plt.plot(x,y2, label='Pol')
plt.plot(x,y3, label='Func')
plt.plot(x,y4, label='N')
plt.plot(x,y5, label='Eps')
plt.xlabel("x")
plt.ylabel("y")
plt.legend()
plt.grid()
plt.show()
