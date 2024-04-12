p=257

x=2
y=2

ans=0
d=x-y

ans+=(p+1)//2*(d**(p-1))%257
print("ans1: ",ans)

for i in range(1,p-1,2):
	print("i: ",i)
	c=0
	f=(d**i)%p
	for a in range(1,(p-1)//2+1):
		c+=(a**(p-1-i))%p
	c=c%p
	print("c: ",c)
	print("temp1: ",f)
	h=(f*c)%p
	print("temp: ",h)
	ans+=h
	print("")

print(ans%p)
