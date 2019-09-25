#!/usr/bin/env python
# coding: utf-8

# In[71]:


import numpy as np
from numpy.linalg import norm
f = open('out_data.txt', 'w')


# In[2]:


A = np.loadtxt("A.txt", delimiter=' ', dtype=np.float)
epsilon = 10 **(-5)
U,L = np.zeros(A.shape), np.zeros(A.shape)
print("Matrix A")
print(A)


# In[3]:


def qr(A):
    m, n = A.shape
    Q = np.eye(m)
    for i in range(n - (m == n)):
        H = np.eye(m)
        H[i:, i:] = make_householder(A[i:, i])
        Q = np.dot(Q, H)
        A = np.dot(H, A)
    return Q, A
 
def make_householder(a):
    v = a / (a[0] + np.copysign(np.linalg.norm(a), a[0]))
    v[0] = 1
    H = np.eye(a.shape[0])
    H -= (2 / np.dot(v, v)) * np.dot(v[:, None], v[None, :])
    return H


# In[4]:


Q,R = qr(A)
print("Matrix q")
print(Q)
print("Matrix r")
print(R)


# In[5]:


A_next = A
iteration = 1
vekt = 0
while(True):
    A_prev = A_next
    q,r = qr(A_prev)
    print("\nIteration:", iteration)
    print("\nMatrix Q\n")
    print(q)
    print('\nMatrix R\n')
    print(r)
    A_next = np.dot(r,q)
    if iteration == 1:
        vekt = q
    else:
        vekt = np.dot(vekt,q)
    err = np.abs(np.diag(A_prev)-np.diag(A_next))
    print("\nCriterion of convergence: {}\n".format(err))
    print("\nFrobenious norm: {}\n".format(np.linalg.norm(A_next)))
    
    iteration += 1
    
    x = np.array(np.dot(A,vekt) - np.diag(A_next)*vekt)
    
    print("Vector of disconnection:\n",x)
 
    if(np.all(x< 0.00001)):
        break;
vekt = np.dot(vekt,q)


# In[7]:


print("Eigenvalues: ",np.diag(A_next))
print("\nEigenvectors:\n\n",vekt)


# In[66]:





# In[67]:


def stepenevmin(a):
    lmax, vmax, imax = stepenevmax(a)
    b = a.copy()
    for i in range(a.shape[0]):
        b[i][i] -= lmax
    x = np.array([1.0] * a.shape[0])
    l2 = 1
    eps = 0.00001
    i= 0
    while True:
        i += 1
        l1 = l2
        xk = np.dot(b, x)
        l2 = xk[0] / x[0]
        norma = norm(xk, ord=2)
        x = xk / norma
        if fabs(l1 - l2) <= eps:
            break
    x = -1*x
    lmin = l2 + lmax
    return lmin, x, i

def out_stepenev(a):
    
    lmax, vmax, imax = stepenevmax(a)
    lmin, vmin, imin = stepenevmin(a)
    f.write("\n Power iteration method:\n")
    f.write("Max eigen value: ")
    f.write(str(lmax))
    f.write("\nVector: \n")
    out_vec_in_file(vmax)
    f.write("Steps: ")
    f.write(str(imax))
    f.write("\n\nMin eigen value: ")
    f.write(str(lmin))
    f.write("\nVector: \n")
    out_vec_in_file(vmin)
    f.write("Steps: ")
    f.write(str(imin))
    f.write("\n")
    f.close()


# In[68]:


def out_vec_in_file(v):
    for i in range(len(v)):
        f.write(str(v[i]) + " ")
    f.write("\n")


# In[69]:


def stepenevmax(a):
    x = np.array([1.0]*a.shape[0])
    l2 = 1
    eps = 0.00001
    i= 0
    while True:
        i += 1
        l1 = l2
        xk = np.dot(a, x)
        l2 = xk[0]/x[0]
        norma = norm(xk, ord=2)
        x = xk/norma
        if fabs(l1 - l2) <= eps:
            break
    x = -1 * x
    return l2, x, i


# In[70]:


out_stepenev(A)


# In[72]:


f.close()


# In[ ]:




