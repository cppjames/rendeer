from PIL import Image, ImageDraw
import math
import meshes
from tqdm import tqdm
import settings

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
    
    for v in tqdm(vertices, ascii = True, desc = "Rotating the mesh"):
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
    
def sortFacesByDistance(arr): 
    for i in tqdm(range(1, len(arr)), ascii = True, desc = "Sorting faces"): 
  
        key = arr[i]
  
        j = i-1
        while j >=0 and key[0] < arr[j][0] : 
                arr[j+1] = arr[j]
                j -= 1
        arr[j+1] = key 

shadingMode = settings.shadingMode
    
yRotation = settings.yRotation
xRotation = settings.xRotation
zRotation = settings.zRotation

lightDirectionX = settings.lightDirectionX
lightDirectionY = settings.lightDirectionY
lightDirectionZ = settings.lightDirectionZ

diffuseColor = settings.diffuseColor

lightVector = settings.lightVector
    
meshPosition = settings.meshPosition

renderFaceNormals = settings.renderFaceNormals
renderHardVertexNormals = settings.renderHardVertexNormals
renderSoftVertexNormals = settings.renderSoftVertexNormals
renderMeshCenter = settings.renderMeshCenter
renderFaceCenters = settings.renderFaceCenters
renderLightSource = settings.renderLightSource

initialVertices = []

for vertex in tqdm(meshes.meshVertices, ascii = True, desc = "Loading vertices"):
    initialVertices.append(vertex)

verticesPosition = list(map(lambda x: add(x, meshPosition), initialVertices))
centerPoint = center(verticesPosition)

verticesOnScreen = []
rotatedVerticesXYZ = []

faceNormals = []

cameraPosition = settings.cameraPosition

imagePlaneZ = settings.imagePlaneZ

xRotation = degreesToRadians(xRotation)
yRotation = degreesToRadians(yRotation)
zRotation = degreesToRadians(zRotation)

lightDirectionX = degreesToRadians(lightDirectionX)
lightDirectionY = degreesToRadians(lightDirectionY)
lightDirectionZ = degreesToRadians(lightDirectionZ)

hardVertexNormals = []

rotatedVerticesXYZ = rotatePoints(verticesPosition, centerPoint, yRotation, xRotation, zRotation)

lightVector = rotatePoints(lightVector, [0, 0, 0], lightDirectionY, lightDirectionX, lightDirectionZ)
    
newLightPoint = calculateScreenPoint(lightVector[0])

for v in tqdm(rotatedVerticesXYZ, ascii = True, desc = "Calculating screen points"):
    
    newPoint = calculateScreenPoint(v)
    
    verticesOnScreen.append(newPoint)

newCenterPoint = calculateScreenPoint(centerPoint)

meshFaces = []
mesh3DFaces = []

for face in tqdm(meshes.meshFaces, ascii = True, desc = "Loading faces"):
    meshFaces.append([verticesOnScreen[face[0]], verticesOnScreen[face[1]], verticesOnScreen[face[2]]])
    mesh3DFaces.append([rotatedVerticesXYZ[face[0]], rotatedVerticesXYZ[face[1]], rotatedVerticesXYZ[face[2]]])
    
gridCenterOnScreen = calculateScreenPoint([0, 0, 0])
    
faceNormalsOnScreen = []
faceCenters = []

for face in tqdm(mesh3DFaces, ascii = True, desc = "Calculating face centers"):
    faceCenters.append(center([face[0], face[1], face[2]]))

sortedDistances = []
c = 0
for center in faceCenters:
    sortedDistances.append([distance(center, cameraPosition), c])
    c += 1
    
sortFacesByDistance(sortedDistances)
sortedDistances.reverse()
tempMeshFaces = []
tempMesh3DFaces = []
for d in sortedDistances:
    tempMeshFaces.append(meshFaces[d[1]])
    tempMesh3DFaces.append(mesh3DFaces[d[1]])
meshFaces = tempMeshFaces
mesh3DFaces = tempMesh3DFaces
    
k = 0
    
for face in tqdm(mesh3DFaces, ascii = True, desc = "Calculating face normals"):
    faceNormal = calculateSurfaceNormal(face)
    
    faceNormals.append([faceNormal[0] + faceCenters[k][0], faceNormal[1] + faceCenters[k][1], faceNormal[2] + faceCenters[k][2]])
    
    k += 1
    
if (settings.calculateVertexNormals):
    for vertex in tqdm(rotatedVerticesXYZ, ascii = True, desc = "Calculating vertex normals"):
        vertexNormals = []
        h = 0
        for face in mesh3DFaces:
            if vertex in face:
                vertexNormals.append([faceNormals[h][0] - faceCenters[h][0] + vertex[0], faceNormals[h][1] - faceCenters[h][1] + vertex[1], faceNormals[h][2] - faceCenters[h][2] + vertex[2]])
            h += 1
        hardVertexNormals.append(vertexNormals)
    
if (renderFaceNormals):
    for n in faceNormals:
        newPoint = calculateScreenPoint(n)
        faceNormalsOnScreen.append(newPoint)
    
faceCentersScreen = []
if (renderFaceCenters):
    for n in faceCenters:
        newPoint = calculateScreenPoint(n)
        faceCentersScreen.append(newPoint)
    
verticesNormalsOnScreen = []
if (renderHardVertexNormals):
    for v in hardVertexNormals:
        vertexNormalsOnScreen = []
        for n in v:
            newPoint = calculateScreenPoint(n)
            vertexNormalsOnScreen.append(newPoint)
        verticesNormalsOnScreen.append(vertexNormalsOnScreen)
facesColors = []
j = 0
for normal in tqdm(faceNormals, ascii = True, desc = "Calculating face colors"):
    lightNormalAngle = angleBetweenVectors([normal[0] - faceCenters[j][0], normal[1] - faceCenters[j][1], normal[2] - faceCenters[j][2]], lightVector[0])
    
    if (lightNormalAngle <= 90):
        if (shadingMode == settings.ShadingMode.flatDiffuse):
            facesColors.append((2, 2, 2))
        elif (shadingMode == settings.ShadingMode.unlit):
            facesColors.append((settings.unlitColor[0], settings.unlitColor[1], settings.unlitColor[2]))
    else:
        if (shadingMode == settings.ShadingMode.flatDiffuse):
            lightNormalAngle -= 90
            color = (int(lightNormalAngle * diffuseColor[0] / 90), int(lightNormalAngle * diffuseColor[1] / 90), int(lightNormalAngle * diffuseColor[2] / 90))
            facesColors.append(color)
        elif (shadingMode == settings.ShadingMode.unlit):
            facesColors.append((settings.unlitColor[0], settings.unlitColor[1], settings.unlitColor[2]))
    
    j += 1

    
img = Image.new('RGB', (300, 300), color = 'black')

draw = ImageDraw.Draw(img)
    
i = 0
for face in tqdm(meshFaces, ascii = True, desc = "Rendering faces"):
    if (distance(faceCenters[i], cameraPosition) > distance(faceNormals[i], cameraPosition)) :
        if (shadingMode == settings.ShadingMode.unlit or shadingMode == settings.ShadingMode.flatDiffuse):
            draw.polygon([coord for vertex in face for coord in vertex], fill = facesColors[i])
        elif (shadingMode == settings.ShadingMode.wireframe):
            draw.polygon([coord for vertex in face for coord in vertex], outline = (settings.wireframeColor[0], settings.wireframeColor[1], settings.wireframeColor[2]))
        if (renderFaceNormals):
            draw.line([faceCentersScreen[i][0], faceCentersScreen[i][1], faceNormalsOnScreen[i][0], faceNormalsOnScreen[i][1]], (160, 0, 0))
    i += 1
    
if (renderHardVertexNormals):
    m = 0
    for vertex in verticesNormalsOnScreen:
        for normal in vertex:
            draw.line([verticesOnScreen[m][0], verticesOnScreen[m][1], normal[0], normal[1]], (0, 0, 200))
        m += 1
    
if (renderFaceCenters):
    for point in faceCentersScreen:
        draw.ellipse((point[0] - 1, point[1] - 1, point[0] + 1, point[1] + 1), fill = 'red')
    
if (renderMeshCenter):
    draw.ellipse((newCenterPoint[0] - 1, newCenterPoint[1] - 1, newCenterPoint[0] + 1, newCenterPoint[1] + 1), fill = 'gray')

if (renderLightSource):
    draw.ellipse((gridCenterOnScreen[0] - 2, gridCenterOnScreen[1] - 2, gridCenterOnScreen[0] + 2, gridCenterOnScreen[1] + 2), fill = 'lightgray')
    draw.line([gridCenterOnScreen[0], gridCenterOnScreen[1], newLightPoint[0], newLightPoint[1]], "yellow")
    
img.save('render.png')
img.show()
