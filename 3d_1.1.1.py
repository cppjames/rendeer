from PIL import Image, ImageDraw
import math
import time
from time import sleep

def add(v1, v2):
    return [x + y for (x, y) in zip(v1, v2)]
    
def substract(v1, v2):
    return [x - y for (x, y) in zip(v1, v2)]
    
def calculateScreenPoint(v):
    xv = v[0]
    yv = v[1]
    zv = v[2]
    
    xk = cameraPosition[0]
    yk = cameraPosition[1]
    zk = cameraPosition[2]
    
    zp = imagePlaneZ
    
    xp = xv - (((xv - xk) * (zv - zp)) / (zv - zk))
    yp = yv - (((yv - yk) * (zv - zp)) / (zv - zk))
    
    return [int(xp * 50 + 100), int(yp * 50 + 100)]
    
def calculateSurfaceNormal(triangle):
    U = substract(triangle[1], triangle[0])
    V = substract(triangle[2], triangle[0])
    
    Nx = (U[1] * V[2]) - (U[2] * V[1])
    Ny = (U[2] * V[0]) - (U[0] * V[2])
    Nz = (U[0] * V[1]) - (U[1] * V[0])
    
    return [Nx/20, Ny/20, Nz/20]
    
def distance(a, b):
    return (math.sqrt(math.pow((b[0] - a[0]), 2) + math.pow((b[1] - a[1]), 2) + math.pow((b[2] - a[2]), 2)))

yRotation = 35 #rotation on the Y axis (in degrees)
xRotation = -25 #rotation on the X axis (in degrees)
zRotation = 0 #rotation on the Z axis (in degrees)
    
meshPosition = [-1, 0, 2] #position of the mesh

initialVertices = [[0, 0, 0], [4, 0, 0], [4, 0, 4], [0, 0, 4], [0, 4, 0], [4, 4, 0], [4, 4, 4], [0, 4, 4]] #mesh vertices
verticesPosition = list(map(lambda x: add(x, meshPosition), initialVertices))
centerPoint = [2 + meshPosition[0], 2 + meshPosition[1], 2 + meshPosition[2]] #center point of the mesh

verticesOnScreen = []
rotatedVerticesY = []
rotatedVerticesXY = []
rotatedVerticesXYZ = []

faceNormals = []

cameraPosition = [1, 1, -2] #position of the camera

imagePlaneZ = 0 #position of the plane (on the Z axis)

xRotation = (xRotation * 2 * math.pi) / 360
yRotation = (yRotation * 2 * math.pi) / 360
zRotation = (zRotation * 2 * math.pi) / 360

for v in verticesPosition:
    v[0] -= centerPoint[0]
    v[1] -= centerPoint[1]
    v[2] -= centerPoint[2]
    rotatedx = v[0]
    rotatedy = v[1]
    rotatedz = v[2]
    
    if (yRotation != 0):
        a = v[0]
        b = v[2]
        rotatedx = a * math.cos(yRotation) + b * math.sin(yRotation)
        rotatedz = b * math.cos(yRotation) - a * math.sin(yRotation)
        rotatedy = v[1]
        
    rotatedx += centerPoint[0]
    rotatedy += centerPoint[1]
    rotatedz += centerPoint[2]
        
    rotatedVerticesY.append([rotatedx, rotatedy, rotatedz])
    
for v in rotatedVerticesY:
    v[0] -= centerPoint[0]
    v[1] -= centerPoint[1]
    v[2] -= centerPoint[2]
    rotatedx = v[0]
    rotatedy = v[1]
    rotatedz = v[2]
    
    if (xRotation != 0):
        a = v[1]
        b = v[2]
        rotatedx = v[0]
        rotatedz = b * math.cos(xRotation) - a * math.sin(xRotation)
        rotatedy = a * math.cos(xRotation) + b * math.sin(xRotation)
        
    rotatedx += centerPoint[0]
    rotatedy += centerPoint[1]
    rotatedz += centerPoint[2]
        
    rotatedVerticesXY.append([rotatedx, rotatedy, rotatedz])
    
for v in rotatedVerticesXY:
    v[0] -= centerPoint[0]
    v[1] -= centerPoint[1]
    v[2] -= centerPoint[2]
    rotatedx = v[0]
    rotatedy = v[1]
    rotatedz = v[2]
    
    if (zRotation != 0):
        a = v[0]
        b = v[1]
        rotatedx = a * math.cos(zRotation) + b * math.sin(zRotation)
        rotatedz = v[2]
        rotatedy = b * math.cos(zRotation) - a * math.sin(zRotation)
        
    rotatedx += centerPoint[0]
    rotatedy += centerPoint[1]
    rotatedz += centerPoint[2]
        
    rotatedVerticesXYZ.append([rotatedx, rotatedy, rotatedz])
    

for v in rotatedVerticesXYZ:
    
    newPoint = calculateScreenPoint(v)
    
    verticesOnScreen.append(newPoint)

newCenterPoint = calculateScreenPoint(centerPoint)


meshFaces = [[verticesOnScreen[0], verticesOnScreen[1], verticesOnScreen[2]],
            [verticesOnScreen[0], verticesOnScreen[3], verticesOnScreen[2]], 
            [verticesOnScreen[4], verticesOnScreen[5], verticesOnScreen[6]], 
            [verticesOnScreen[4], verticesOnScreen[7], verticesOnScreen[6]],
            [verticesOnScreen[0], verticesOnScreen[1], verticesOnScreen[5]],
            [verticesOnScreen[0], verticesOnScreen[4], verticesOnScreen[5]],
            [verticesOnScreen[1], verticesOnScreen[2], verticesOnScreen[6]],
            [verticesOnScreen[1], verticesOnScreen[5], verticesOnScreen[6]],
            [verticesOnScreen[2], verticesOnScreen[6], verticesOnScreen[7]],
            [verticesOnScreen[2], verticesOnScreen[3], verticesOnScreen[7]],
            [verticesOnScreen[0], verticesOnScreen[4], verticesOnScreen[3]],
            [verticesOnScreen[4], verticesOnScreen[7], verticesOnScreen[3]]]
#lines that connect vertices
#meshLines = [[verticesOnScreen[0], verticesOnScreen[1]], [verticesOnScreen[1], verticesOnScreen[2]], [verticesOnScreen[2], verticesOnScreen[3]], [verticesOnScreen[3], verticesOnScreen[0]], [verticesOnScreen[4], verticesOnScreen[5]], [verticesOnScreen[5], verticesOnScreen[6]], [verticesOnScreen[6], verticesOnScreen[7]], [verticesOnScreen[4], verticesOnScreen[7]], [verticesOnScreen[0], verticesOnScreen[4]], [verticesOnScreen[1], verticesOnScreen[5]], [verticesOnScreen[2], verticesOnScreen[6]], [verticesOnScreen[3], verticesOnScreen[7]]]
mesh3DFaces = [[rotatedVerticesXYZ[0], rotatedVerticesXYZ[1], rotatedVerticesXYZ[2]],
            [rotatedVerticesXYZ[0], rotatedVerticesXYZ[3], rotatedVerticesXYZ[2]], 
            [rotatedVerticesXYZ[4], rotatedVerticesXYZ[5], rotatedVerticesXYZ[6]], 
            [rotatedVerticesXYZ[4], rotatedVerticesXYZ[7], rotatedVerticesXYZ[6]],
            [rotatedVerticesXYZ[0], rotatedVerticesXYZ[1], rotatedVerticesXYZ[5]],
            [rotatedVerticesXYZ[0], rotatedVerticesXYZ[4], rotatedVerticesXYZ[5]],
            [rotatedVerticesXYZ[1], rotatedVerticesXYZ[2], rotatedVerticesXYZ[6]],
            [rotatedVerticesXYZ[1], rotatedVerticesXYZ[5], rotatedVerticesXYZ[6]],
            [rotatedVerticesXYZ[2], rotatedVerticesXYZ[6], rotatedVerticesXYZ[7]],
            [rotatedVerticesXYZ[2], rotatedVerticesXYZ[3], rotatedVerticesXYZ[7]],
            [rotatedVerticesXYZ[0], rotatedVerticesXYZ[4], rotatedVerticesXYZ[3]],
            [rotatedVerticesXYZ[4], rotatedVerticesXYZ[7], rotatedVerticesXYZ[3]]]
            
faceNormalsInvert = [0, 1, 1, 0, 1, 0, 1, 0, 0, 1, 1, 1]

img = Image.new('RGB', (300, 300), color = 'black')

draw = ImageDraw.Draw(img)
#for point in verticesOnScreen:
    #draw.ellipse((point[0] - 2, point[1] - 2, point[0] + 2, point[1] + 2), fill = 'lightgray')
    
    
#for line in meshLines:
    #draw.line((line[0][0], line[0][1], line[1][0], line[1][1]), fill = 'lightgray', width = 1)
    
faceNormalsOnScreen = []
faceCenters = []

for face in mesh3DFaces:
    faceCenters.append([(face[0][0] + face[1][0] + face[2][0])/3, (face[0][1] + face[1][1] + face[2][1])/3, (face[0][2] + face[1][2] + face[2][2])/3])

k = 0
    
for face in mesh3DFaces:
    faceNormal = calculateSurfaceNormal(face)
    
    faceNormals.append([faceNormal[0] + faceCenters[k][0], faceNormal[1] + faceCenters[k][1], faceNormal[2] + faceCenters[k][2]])
    
    k += 1
for n in faceNormals:
    
    newPoint = calculateScreenPoint(n)
    faceNormalsOnScreen.append(newPoint)
    
faceCentersScreen = []
for n in faceCenters:
    
    newPoint = calculateScreenPoint(n)
    faceCentersScreen.append(newPoint)
  
    
i = 0

for face in meshFaces:
    if (distance(faceCenters[i], cameraPosition) > distance(faceNormals[i], cameraPosition)) :
        if (faceNormalsInvert[i] == 0):
            draw.polygon([coord for vertex in face for coord in vertex], outline = "white")
            #draw.line([faceCentersScreen[i][0], faceCentersScreen[i][1], faceNormalsOnScreen[i][0], faceNormalsOnScreen[i][1]], "red")
    else :
        if (faceNormalsInvert[i] == 1):
            draw.polygon([coord for vertex in face for coord in vertex], outline = "white")
    i += 1
    
#for point in faceNormalsOnScreen:
    #draw.ellipse((point[0] - 1, point[1] - 1, point[0] + 1, point[1] + 1), fill = 'red')
    

j = 0  
for point in faceCentersScreen:
    #draw.ellipse((point[0] - 1, point[1] - 1, point[0] + 1, point[1] + 1), fill = 'red')
    
    j += 1
    
draw.ellipse((newCenterPoint[0] - 1, newCenterPoint[1] - 1, newCenterPoint[0] + 1, newCenterPoint[1] + 1), fill = 'gray')
    
img.save('render.png')
img.show()
