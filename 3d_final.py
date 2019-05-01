from PIL import Image, ImageDraw
import math

def add(v1, v2):
    return [x + y for (x, y) in zip(v1, v2)]
    
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
    
    return [xp, yp]

yRotation = 30 #rotation on the Y axis (in degrees)
xRotation = 0 #rotation on the X axis (in degrees)
zRotation = 0 #rotation on the Z axis (in degrees)
    
meshPosition = [-1, 0, 2] #position of the mesh

initialVertices = [[0, 0, 0], [4, 0, 0], [4, 0, 4], [0, 0, 4], [0, 4, 0], [4, 4, 0], [4, 4, 4], [0, 4, 4]] #mesh vertices
verticesPosition = list(map(lambda x: add(x, meshPosition), initialVertices))
centerPoint = [2 + meshPosition[0], 2 + meshPosition[1], 2 + meshPosition[2]] #center point of the mesh

verticesOnScreen = []
rotatedVerticesY = []
rotatedVerticesXY = []
rotatedVerticesXYZ = []

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
    
    verticesOnScreen.append([int(newPoint[0] * 50 + 100), int(newPoint[1] * 50 + 100)])

newCenterPoint = [int(calculateScreenPoint(centerPoint)[0] * 50 + 100), int(calculateScreenPoint(centerPoint)[1] * 50 + 100)]

#lines that connect vertices
meshLines = [[verticesOnScreen[0], verticesOnScreen[1]], [verticesOnScreen[1], verticesOnScreen[2]], [verticesOnScreen[2], verticesOnScreen[3]], [verticesOnScreen[3], verticesOnScreen[0]], [verticesOnScreen[4], verticesOnScreen[5]], [verticesOnScreen[5], verticesOnScreen[6]], [verticesOnScreen[6], verticesOnScreen[7]], [verticesOnScreen[4], verticesOnScreen[7]], [verticesOnScreen[0], verticesOnScreen[4]], [verticesOnScreen[1], verticesOnScreen[5]], [verticesOnScreen[2], verticesOnScreen[6]], [verticesOnScreen[3], verticesOnScreen[7]]]

img = Image.new('RGB', (300, 300), color = 'black')

draw = ImageDraw.Draw(img)
for point in verticesOnScreen:
    draw.ellipse((point[0] - 2, point[1] - 2, point[0] + 2, point[1] + 2), fill = 'lightgray')
    
for line in meshLines:
    draw.line((line[0][0], line[0][1], line[1][0], line[1][1]), fill = 'lightgray', width = 1)
    
draw.ellipse((newCenterPoint[0] - 1, newCenterPoint[1] - 1, newCenterPoint[0] + 1, newCenterPoint[1] + 1), fill = 'gray')
    
img.save('render.png')
img.show()
