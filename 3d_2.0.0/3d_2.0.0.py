from PIL import Image, ImageDraw
import math
import meshes

def center(vertices):
    centerx = 0
    centery = 0
    centery = 0
    sumx = 0
    sumy = 0
    sumz = 0
    for v in vertices:
        sumx += v[0]
        sumy += v[1]
        sumz += v[2]
    centerx = sumx / len(vertices)
    centery = sumy / len(vertices)
    centerz = sumz / len(vertices)
        
    return [centerx, centery, centerz]
        
def rotatePoints(vertices, centerVector, yRot, xRot, zRot):
    rotVerticesY = []
    rotVerticesXY = []
    rotVerticesXYZ = []
    
    for v in vertices:
        v[0] -= centerVector[0]
        v[1] -= centerVector[1]
        v[2] -= centerVector[2]
        rotatedx = v[0]
        rotatedy = v[1]
        rotatedz = v[2]
        
        if (yRot != 0):
            a = v[0]
            b = v[2]
            rotatedx = a * math.cos(yRot) + b * math.sin(yRot)
            rotatedz = b * math.cos(yRot) - a * math.sin(yRot)
            rotatedy = v[1]
            
        rotatedx += centerVector[0]
        rotatedy += centerVector[1]
        rotatedz += centerVector[2]
            
        rotVerticesY.append([rotatedx, rotatedy, rotatedz])
        
    for v in rotVerticesY:
        v[0] -= centerVector[0]
        v[1] -= centerVector[1]
        v[2] -= centerVector[2]
        rotatedx = v[0]
        rotatedy = v[1]
        rotatedz = v[2]
        
        if (xRot != 0):
            a = v[1]
            b = v[2]
            rotatedx = v[0]
            rotatedz = b * math.cos(xRot) - a * math.sin(xRot)
            rotatedy = a * math.cos(xRot) + b * math.sin(xRot)
            
        rotatedx += centerVector[0]
        rotatedy += centerVector[1]
        rotatedz += centerVector[2]
            
        rotVerticesXY.append([rotatedx, rotatedy, rotatedz])
        
    for v in rotVerticesXY:
        v[0] -= centerVector[0]
        v[1] -= centerVector[1]
        v[2] -= centerVector[2]
        rotatedx = v[0]
        rotatedy = v[1]
        rotatedz = v[2]
        
        if (zRot != 0):
            a = v[0]
            b = v[1]
            rotatedx = a * math.cos(zRot) + b * math.sin(zRot)
            rotatedz = v[2]
            rotatedy = b * math.cos(zRot) - a * math.sin(zRot)
            
        rotatedx += centerVector[0]
        rotatedy += centerVector[1]
        rotatedz += centerVector[2]
            
        rotVerticesXYZ.append([rotatedx, rotatedy, rotatedz])
        
    return rotVerticesXYZ

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
    
    return [Nx/5, Ny/5, Nz/5]
    
def distance(a, b):
    return (math.sqrt(math.pow((b[0] - a[0]), 2) + math.pow((b[1] - a[1]), 2) + math.pow((b[2] - a[2]), 2)))
    
def degreesToRadians(angle):
    return (angle * 2 * math.pi) / 360
    
def radiansToDegrees(angle):
    return (angle * 360) / (math.pi * 2)

def angleBetweenVectors(a, b):
    cos = (a[0] * b[0] + a[1] * b[1] + a[2] * b[2]) / (distance([0, 0, 0], a) * distance([0, 0, 0], b))
    return radiansToDegrees(math.acos(cos))
    
yRotation = 0 #rotation on the Y axis (in degrees)
xRotation = 0 #rotation on the X axis (in degrees)
zRotation = -10 #rotation on the Z axis (in degrees)

lightDirectionX = 0
lightDirectionY = -50
lightDirectionZ = -20

lightVector = [[1, 0, 0]]
randomVector = [0, 0, 1]
    
meshPosition = [1.5, 1, 0] #position of the mesh

initialVertices = []

for vertex in meshes.sphereVertices:
    initialVertices.append(vertex)

verticesPosition = list(map(lambda x: add(x, meshPosition), initialVertices))
centerPoint = center(verticesPosition) #center point of the mesh

verticesOnScreen = []
rotatedVerticesXYZ = []

faceNormals = []

cameraPosition = [1, 1, -2] #position of the camera

imagePlaneZ = 0 #position of the plane (on the Z axis)

xRotation = degreesToRadians(xRotation)
yRotation = degreesToRadians(yRotation)
zRotation = degreesToRadians(zRotation)

lightDirectionX = degreesToRadians(lightDirectionX)
lightDirectionY = degreesToRadians(lightDirectionY)
lightDirectionZ = degreesToRadians(lightDirectionZ)

rotatedVerticesXYZ = rotatePoints(verticesPosition, centerPoint, yRotation, xRotation, zRotation)

lightVector = rotatePoints(lightVector, [0, 0, 0], lightDirectionY, lightDirectionX, lightDirectionZ)
    
newLightPoint = calculateScreenPoint(lightVector[0])

for v in rotatedVerticesXYZ:
    
    newPoint = calculateScreenPoint(v)
    
    verticesOnScreen.append(newPoint)

newCenterPoint = calculateScreenPoint(centerPoint)

meshFaces = []
mesh3DFaces = []

for face in meshes.sphereFaces:
    meshFaces.append([verticesOnScreen[face[0]], verticesOnScreen[face[1]], verticesOnScreen[face[2]]])
    mesh3DFaces.append([rotatedVerticesXYZ[face[0]], rotatedVerticesXYZ[face[1]], rotatedVerticesXYZ[face[2]]])
    
faceNormalsInvert = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]

img = Image.new('RGB', (300, 300), color = 'black')

draw = ImageDraw.Draw(img)
    
gridCenterOnScreen = calculateScreenPoint([0, 0, 0])
    
faceNormalsOnScreen = []
faceCenters = []

for face in mesh3DFaces:
    faceCenters.append(center([face[0], face[1], face[2]]))

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
  

facesColors = []
j = 0
for normal in faceNormals:
    lightNormalAngle = angleBetweenVectors([normal[0] - faceCenters[j][0], normal[1] - faceCenters[j][1], normal[2] - faceCenters[j][2]], lightVector[0])
    
    if (lightNormalAngle <= 90):
        facesColors.append((40, 40, 40))
    else:
        lightNormalAngle -= 90
        color = int(lightNormalAngle * 255 / 90)
        facesColors.append((color, color, color))
    
    j += 1

i = 0
for face in meshFaces:
    if (distance(faceCenters[i], cameraPosition) > distance(faceNormals[i], cameraPosition)) :
        if (faceNormalsInvert[i] == 0):
            draw.polygon([coord for vertex in face for coord in vertex], fill = facesColors[i])
            draw.line([faceCentersScreen[i][0], faceCentersScreen[i][1], faceNormalsOnScreen[i][0], faceNormalsOnScreen[i][1]], (160, 0, 0))
    else :
        if (faceNormalsInvert[i] == 1):
            draw.polygon([coord for vertex in face for coord in vertex], outline = "white", fill = "gray")
    i += 1
    

#for point in faceCentersScreen:
    #draw.ellipse((point[0] - 1, point[1] - 1, point[0] + 1, point[1] + 1), fill = 'red')
    
draw.ellipse((newCenterPoint[0] - 1, newCenterPoint[1] - 1, newCenterPoint[0] + 1, newCenterPoint[1] + 1), fill = 'gray')
#draw.ellipse((gridCenterOnScreen[0] - 2, gridCenterOnScreen[1] - 2, gridCenterOnScreen[0] + 2, gridCenterOnScreen[1] + 2), fill = 'lightgray')
#draw.line([gridCenterOnScreen[0], gridCenterOnScreen[1], newLightPoint[0], newLightPoint[1]], "yellow")
    
img.save('render.png')
img.show()
