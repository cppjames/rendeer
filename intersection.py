import math

S = []
P = [None] * 6
I = []

P[0] = [60, 400];
P[1] = [150, 80];
P[2] = [500, 400];
P[3] = [700, 200];

P[4] = [300, 50];
P[5] = [400, 450];


def sgn( x ):
    if (x < 0.0):
        return -1
    else:
        return 1

def abs(x):
    if (x < 0.0):
        return -x
    else:
        return x
        
def cubicRoots(P):
    a=P[0]
    b=P[1]
    c=P[2]
    d=P[3]
    
    A=b/a
    B=c/a
    C=d/a
    
    Q = []
    R = []
    D = []
    S = []
    T = []
    Im = []

    Q = (3*B - A**2)/9
    R = (9*A*B - 27*C - 2*A**3)/54
    D = Q**3 + R**2

    t = [None] * 3
	
    if (D >= 0):
        S = sgn(R + D**0.5)*abs(R + D**0.5)**(1/3)
        T = sgn(R - D**0.5)*abs(R - D**0.5)**(1/3)

        t[0] = -A/3 + (S + T)
        t[1] = -A/3 - (S + T)/2
        t[2] = -A/3 - (S + T)/2
        Im = abs((3**0.5)*(S - T)/2)
        
        if (Im!=0):
            t[1]=-1
            t[2]=-1
    
    else:
        th = math.acos(R/(-Q**3)**0.5)
        
        t[0] = 2*math.sqrt(-Q)*math.cos(th/3) - A/3
        t[1] = 2*math.sqrt(-Q)*math.cos((th + 2*math.pi)/3) - A/3
        t[2] = 2*math.sqrt(-Q)*math.cos((th + 4*math.pi)/3) - A/3
        Im = 0.0
    
    for i in range(3):
        if (t[i]<0 or t[i]>1.0):
            t[i]=-1
                
    t=sortSpecial(t)
    
    return t

def sortSpecial(a):
    flip=False;
    for i in range(len(a)-1):
        if ((a[i+1]>=0 and a[i]>a[i+1]) or (a[i]<0 and a[i+1]>=0)):
            flip=True
            temp=a[i]
            a[i]=a[i+1]
            a[i+1]=temp
    while (flip):
        flip=False;
        for i in range(len(a)-1):
            if ((a[i+1]>=0 and a[i]>a[i+1]) or (a[i]<0 and a[i+1]>=0)):
                flip=True
                temp=a[i]
                a[i]=a[i+1]
                a[i+1]=temp
    return a

def bezierCoeffs(P0,P1,P2,P3):
    Z = [None] * 4
    Z[0] = -P0 + 3*P1 + -3*P2 + P3
    Z[1] = 3*P0 - 6*P1 + 3*P2
    Z[2] = -3*P0 + 3*P1
    Z[3] = P0
    return Z

def computeIntersections(px,py,lx,ly):
    X = [None] * 2;
    point = []
    A=ly[1]-ly[0]
    B=lx[0]-lx[1]
    C=lx[0]*(ly[0]-ly[1]) + ly[0]*(lx[1]-lx[0])
    
    bx = bezierCoeffs(px[0],px[1],px[2],px[3])
    by = bezierCoeffs(py[0],py[1],py[2],py[3])
	
    P = [None] * 4
    P[0] = A*bx[0]+B*by[0]
    P[1] = A*bx[1]+B*by[1]
    P[2] = A*bx[2]+B*by[2]
    P[3] = A*bx[3]+B*by[3] + C
    
    r=cubicRoots(P);
    
    
    for i in range(3):
        t=r[i];
        
        X[0]=bx[0]*t*t*t+bx[1]*t*t+bx[2]*t+bx[3];
        X[1]=by[0]*t*t*t+by[1]*t*t+by[2]*t+by[3];            
            
        if ((lx[1]-lx[0])!=0):
            s=(X[0]-lx[0])/(lx[1]-lx[0])
        else:
            s=(X[1]-ly[0])/(ly[1]-ly[0])
        return([X])
    return point