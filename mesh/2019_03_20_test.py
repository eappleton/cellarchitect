import numpy as np

point1=[5,0,0]
point2=[0,5,0]
point3=[1,1,5]
point4=[-1,-1,0]
cross=np.cross(np.array([point2[0]-point1[0], point2[1]-point1[1], point2[2]-point1[2]]),
         np.array([point3[0]-point1[0], point3[1]-point1[1], point3[2]-point1[2]]))
factor=((point1[0]-point4[0])*cross[0]+(point1[1]-point4[1])*cross[1]+(point1[2]-point4[2])*cross[2])/(cross[0]**2+cross[1]**2+cross[2]**2)
if factor<0:
    temp=point2
    point2=point3
    point3=temp
print(point1,point2,point3,point4)