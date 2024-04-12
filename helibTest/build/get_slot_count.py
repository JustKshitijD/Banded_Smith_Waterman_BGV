p=43
m=490025
slot_count=0

def gcd(a, b):
 
    if (a == 0):
        return b
    return gcd(b % a, a)
 
# A simple method to evaluate
# Euler Totient Function
def phi(n):
 
    result = 1
    for i in range(2, n):
        if (gcd(i, n) == 1):
            result+=1
    return result


def get_ord_p(p,m):
	k=1
	cnt=1
	while (k*p)%m!=1:
		cnt+=1
		k*=p
	return cnt

while True:
	# ord(p) is minimum k such that p^k=1(mod m)
	ord_p=get_ord_p(p,m)
	phi_m=phi(m)
	slot_count=phi_m//ord_p
	print("m: ",m)
	print("slot_count: ",slot_count)
	print(" ")
	if slot_count>=1000:
		break
	m+=1

print("m: ",m)
print("slot_count: ",slot_count)

# m=150444; slot_count: 2376
