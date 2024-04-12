qlen=[1,2,3,4]
max_s=[1,2,3,4,5]

end_bonus=5
o_ins=8
e_ins=3
w=2

interpolation_points=[]

for i in range(0,len(qlen)):
	for j in range(0,len(max_s)):
		x=qlen[i]*max_s[j]+end_bonus-o_ins
		y=min(w,max(1,(int)(x*1.0/3+1)))
		interpolation_points.append((x,y))

print(interpolation_points)

