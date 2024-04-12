import sympy
from sympy import symbols, I, expand
import numpy as np
import sys

from sympy import symbols, Poly


def get_modular_inverse(x):
	global p
	for i in range(1,p):
		if (x*i)%p==1:
			# print("Modular inverse of x="+str(x)+" modulo "+str(p)+" is: "+str(i))
			break
	return i


x = sympy.Symbol('x')

# print((((x-2)/-1*(x-3)/-2*(x-4)/-3)+(x-1)*(x-3)/-1*(x-4)*-1+(x-1)/2*(x-2)*(x-4)*-2+(x-1)/3*(x-2)/2*(x-3)*2).expand())

#######################

# qlen in {1,2,3,4}
# max_s in Enc{1,2,3,4,5}
# eb=5, o_ins=8, e_ins=3
# val = min(w,max(1,(qlen*max_s+eb-o_ins)/e_ins+1))

# ((qlen*max_s+eb-o_ins),val) = [(-2, 1), (-1, 1), (0, 1), (1, 1), (2, 1), (-1, 1), 
# (1, 1), (3, 2), (5, 2), (7, 2), (0, 1), (3, 2), (6, 2), (9, 2), (12, 2), (1, 1), 
# (5, 2), (9, 2), (13, 2), (17, 2)]

#######################

p=43
# list of points to be interpolated
# l= [(-2, 1), (-1, 1), (0, 1), (1, 1), (2, 1), (3, 2), (5, 2), (7, 2), (6, 2), (9, 2), (12, 2), (13, 2), (17, 2)]

qlen=[1,2,3,4]
max_s=[1,2,3,4,5]		# encrypted
eb=1
o_ins=0
e_ins=2
w=4
l=[]
s=set()

for i in range(0,len(qlen)):
	for j in range(0,len(max_s)):
		s.add((qlen[i]*max_s[j]+eb-o_ins))

s=list(s)
print("s: ",s)

for k in s:
	y=min(w,max(1,(int)(k*1.0/e_ins+1)))
	l.append((k,y))

print("l: ",l)

res=0
dens=[]

# get all denominators in the lagrange polynomial expression
for i in range(0,len(l)):
	z=1
	for j in range(0,len(l)):
		if i!=j:
			z*=((l[i][0]-l[j][0]))
	dens.append(z)

# find modular inverse of each denominator
inv_dens= list(map(get_modular_inverse ,dens))

# get all numerators of lagrange expression and represent in single sympy expression
for i in range(0,len(l)):
	z=1
	for j in range(0,len(l)):
		if i!=j:
			z*=((x-l[j][0]))
	z*=(l[i][1]*inv_dens[i])%p
	res+=z	

print(res.expand())

# Get the coefficients of the polynomial
coefficients = Poly(res, x).all_coeffs()

# Find the modulus of each coefficient mod p
modulus_coefficients = [coef%p for coef in coefficients]

# Print the result
print("Original coefficients:", coefficients)
print("Modulus coefficients:", modulus_coefficients)

# Verify the modulus coefficient expression got
# Here, if modulus coefficient=[2,4,5], then, the lagrnage interpolated polynomial is 2x^2+4x+5
# and we evaluate it at x=2
print("Result... - ",np.polyval(modulus_coefficients,2)%p)
# Verify the sympy polynomial before modulus of each coefficient was taken
print(res.evalf(subs={x:2})%p)
