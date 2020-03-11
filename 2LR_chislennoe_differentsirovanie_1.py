import numpy as np 
import matplotlib.pyplot as plt 

def f(x):
    return np.log(np.cosh(x))

def df(x):
    return np.sinh(x)/np.cosh(x)

def d2f(x):
    return -np.sinh(x)*np.sinh(x)/np.cosh(x)/np.cosh(x)+1
def d3f(x):
    return 2*np.sinh(x)/np.cosh(x)*(np.sinh(x)*np.sinh(x)/np.cosh(x)/np.cosh(x)-1)

h = 0.01
a = -3
b = 3
N = int((b-a)/h+1)
X = np.linspace(a ,b, N)
Y = f(X)
dY = df(X)
# Правая разность
rddx = np.zeros(N)
i = 0
while i < N-1:
    rddx[i] = (Y[i+1]-Y[i])/h
    i += 1
rddx[N-1] = (Y[N-1]-Y[N-2])/h

fig = plt.figure(figsize = (15,15))
ax = fig.add_subplot(2,2,1)
ax.grid()
plt.title('Правая разность для шага h = {}'.format(h))
plt.xlabel('x')
plt.ylabel('dy/dx')
plt.plot(X, rddx, 'g')
#plt.plot(X, dY, 'b')
#Центральная разность
mddx = np.zeros(N)
i = 1
while i < N-1:
    mddx[i] = 0.5*(Y[i+1]-Y[i-1])/h
    i += 1
mddx[0] = 0.5*(-3*Y[0]+4*Y[1]-Y[2])/h
mddx[N-1] = 0.5*(3*Y[N-1]-4*Y[N-2]+Y[N-3])/h

ax = fig.add_subplot(2,2,3)
ax.grid()
plt.title('Центральная разность для шага h = {}'.format(h))
plt.xlabel('x')
plt.ylabel('dy/dx')
plt.plot(X, mddx, 'g')
#plt.plot(X, dY, 'b')

h = pow(10,-4)
E = []
H = []
Em = []
while h < 1:
    i = 0
    N = int((b-a)/h+1)
    X = np.linspace(a ,b, N)
    rddx = np.zeros(N)
    mddx = np.zeros(N)
    error = np.zeros(N)
    merror = np.zeros(N)
    dY = df(X)
    Y = f(X)
    while i < N-1:
        rddx[i] = ((Y[i+1])-Y[i])/h
        i += 1
    rddx[N-1] = (Y[N-1]-Y[N-2])/h
    i = 0
    while i < N:
        error[i] = np.abs(rddx[i]-dY[i])
        i += 1
    i = 1
    while i < N-1:
        mddx[i] = 0.5*(Y[i+1]-Y[i-1])/h
        i += 1
    mddx[0] = 0.5*(-3*Y[0]+4*Y[1]-Y[2])/h
    mddx[N-1] = 0.5*(3*Y[N-1]-4*Y[N-2]+Y[N-3])/h
    i = 0
    while i < N:
        merror[i] = np.abs(mddx[i]-dY[i])
        i += 1
    E.append(np.log(max(error)))
    Em.append(np.log(max(merror)))
    H.append(np.log(h))
    h *= 10

ax = fig.add_subplot(2,2,2)
ax.grid()
plt.title('Невязка в ln масштабе O(h)')
plt.xlabel('ln(h)')
plt.ylabel('ln(maxerror)')
plt.plot(H, E, 'r')

ax = fig.add_subplot(2,2,4)
ax.grid()
plt.title('Невязка в ln масштабе O(h2)')
plt.xlabel('ln(h)')
plt.ylabel('ln(maxerror)')
plt.plot(H, Em, 'r')

#Вторая производная, O(h2)
h = 0.01
a = -3
b = 3
N = int((b-a)/h+1)
X = np.linspace(a ,b, N)
Y = f(X)
dY = df(X)
d2Y = d2f(X)
d2dx2 = np.zeros(N)
i = 1
while i < N-1:
    d2dx2[i] = (Y[i+1]-2*Y[i]+Y[i-1])/h/h
    i += 1
d2dx2[0] = (2*Y[0]-5*Y[1]+4*Y[2]-Y[3])/h/h
d2dx2[N-1] = (-Y[N-4]+4*Y[N-3]-5*Y[N-2]+2*Y[N-1])/h/h
Er = np.zeros(N)
Er = np.abs(d2dx2-d2Y)

fig = plt.figure(figsize = (15,15))
ax = fig.add_subplot(3,2,1)
ax.grid()
plt.title('Погрешность для шага h = {} О(h2)'.format(h))
plt.xlabel('x')
plt.ylabel('Error')
plt.plot(X, Er, 'r')

d2dx4 = np.zeros(N)
i = 2
while i < N-2:
    d2dx4[i] = (-Y[i+2]+16*Y[i+1]-30*Y[i]+16*Y[i-1]-Y[i-2])/h/h/12
    i += 1
d2dx4[0] = (-Y[2]+16*Y[1]-30*Y[0]+16*f(X[0]-h)-f(X[0]-2*h))/h/h/12
d2dx4[1] = (-Y[3]+16*Y[2]-30*Y[1]+16*Y[0]-f(X[0]-h))/h/h/12
d2dx4[N-2] = (-f(X[N-1]+h)+16*Y[N-1]-30*Y[N-2]+16*Y[N-3]-Y[N-4])/h/h/12
d2dx4[N-1] = (-f(X[N-1]+2*h)+16*f(X[N-1]+h)-30*Y[N-1]+16*Y[N-2]-Y[N-3])/h/h/12
Er2 = np.zeros(N)
Er2 = np.abs(d2dx4-d2Y)

ax = fig.add_subplot(3,2,3)
ax.grid()
plt.title('Погрешность для шага h = {} О(h4)'.format(h))
plt.xlabel('x')
plt.ylabel('Error')
plt.plot(X, Er2, 'r')


#Ln невязка для второй производной
h = pow(10,-3)
E = []
H = []
Em = []
while h < 1:
    N = int((b-a)/h+1)
    X = np.linspace(a ,b, N)
    Y = f(X)
    dY = df(X)
    d2Y = d2f(X)
    d2dx2 = np.zeros(N)
    d2dx4 = np.zeros(N)
    error2 = np.zeros(N)
    error4 = np.zeros(N)
    i = 1
    while i < N-1:
        d2dx2[i] = (Y[i+1]-2*Y[i]+Y[i-1])/h/h
        i += 1
    d2dx2[0] = (2*Y[0]-5*Y[1]+4*Y[2]-Y[3])/h/h
    d2dx2[N-1] = (-Y[N-4]+4*Y[N-3]-5*Y[N-2]+2*Y[N-1])/h/h
    
    i = 0
    while i < N:
        error2[i] = np.abs(d2dx2[i]-d2Y[i])
        i += 1
    
    i = 2
    while i < N-2:
        d2dx4[i] = (-Y[i+2]+16*Y[i+1]-30*Y[i]+16*Y[i-1]-Y[i-2])/h/h/12
        i += 1
    d2dx4[0] = (-Y[2]+16*Y[1]-30*Y[0]+16*f(X[0]-h)-f(X[0]-2*h))/h/h/12
    d2dx4[1] = (-Y[3]+16*Y[2]-30*Y[1]+16*Y[0]-f(X[0]-h))/h/h/12
    d2dx4[N-2] = (-f(X[N-1]+h)+16*Y[N-1]-30*Y[N-2]+16*Y[N-3]-Y[N-4])/h/h/12
    d2dx4[N-1] = (-f(X[N-1]+2*h)+16*f(X[N-1]+h)-30*Y[N-1]+16*Y[N-2]-Y[N-3])/h/h/12
    
    i = 0
    while i < N:
        error4[i] = np.abs(d2dx4[i]-d2Y[i])
        i += 1
    E.append(np.log(max(error2)))
    Em.append(np.log(max(error4)))
    H.append(np.log(h))
    h *= 10

ax = fig.add_subplot(3,2,2)
ax.grid()
plt.title('ln невязка 2 порядок')
plt.xlabel('ln(h)')
plt.ylabel('Ln(maxer)')
plt.plot(H, E, 'r')

ax = fig.add_subplot(3,2,4)
ax.grid()
plt.title('ln невязка 4 порядок')
plt.xlabel('ln(h)')
plt.ylabel('Ln(maxer)')
plt.plot(H, Em, 'r')

#Третья производная, второй порядок
h = pow(10,-3)
E = []
H = []
while h < 1:
    N = int((b-a)/h+1)
    X = np.linspace(a ,b, N)
    Y = f(X)
    d3Y = d3f(X)
    d3dx2 = np.zeros(N)
    error2 = np.zeros(N)
    i = 2
    while i < N-2:
        d3dx2[i] = 0.5*(Y[i+2]-Y[i-2]+2*Y[i-1]-2*Y[i+1])/h/h/h
        i += 1
    d3dx2[0] = 0.5*(Y[2]-f(X[0]-2*h)+2*f(X[0]-h)-2*Y[1])/h/h/h
    d3dx2[1] = 0.5*(Y[3]-f(X[0]-h)+2*Y[0]-2*Y[2])/h/h/h
    d3dx2[N-2] = 0.5*(f(X[N-1]+h)-Y[N-4]+2*Y[N-3]-2*Y[N-1])/h/h/h
    d3dx2[N-1] = 0.5*(f(X[N-1]+2*h)-Y[N-3]+2*Y[N-2]-2*f(X[N-1]+h))/h/h/h
    i = 0
    while i < N:
        error2[i] = np.abs(d3dx2[i]-d3Y[i])
        i += 1

    E.append(np.log(max(error2)))
    H.append(np.log(h))
    h *= 10

ax = fig.add_subplot(3,2,5)
ax.grid()
plt.title('3 производная ln невязка 2 порядок')
plt.xlabel('ln(h)')
plt.ylabel('Ln(maxer)')
plt.plot(H, E, 'r')