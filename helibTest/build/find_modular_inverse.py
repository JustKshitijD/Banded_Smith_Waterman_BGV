import sys

p=43

x=(int)(sys.argv[1])
print(x)


y=-1

for i in range(1,p):
	if (x*i)%p==1:
		print("Modular inverse of x="+str(x)+" modulo "+str(p)+" is: "+str(i))
		break


# lst = [1126343522304000, 927031705600, 1126343522304000, 75089568153600, 375447840768000, 25029856051200, 86641809408000, 25029856051200, 33127750656000, 243797299200, 15643660032000, 37246809600 , 1]

# fin=[]

# for j in range(0,len(lst)):
# 	x=lst[j]
# 	for i in range(1,p):
# 		if (x*i)%p==1:
# 			# print("Modular inverse of x="+str(x)+" modulo "+str(p)+" is: "+str(i))
# 			fin.append(i)
# 			break
# print(fin)
