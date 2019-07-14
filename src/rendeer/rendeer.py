# -*- coding: utf-8 -*-


import math

from PIL import Image, ImageDraw
from tqdm import tqdm

import rendeer.meshes as meshes
import rendeer.settings as settings
import rendeer.intersection as intersection


# argparse


def main(args=None):
    print(f"Hello from {__name__}")


def center(vertices):
    """Returns the center of a polygon using the arithmetic mean
    """
    sumx, sumy, sumz = 0, 0, 0
    length = len(vertices)
    for v in vertices:
        sumx += v[0]
        sumy += v[1]
        sumz += v[2]
    centerx = sumx / length
    centery = sumy / length
    centerz = sumz / length

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

        if yRot != 0:
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

        if xRot != 0:
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

        if zRot != 0:
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

    xp = xv - (((xv - xk) * (zv - zp)) / ((zv - zk) if (zv - zk != 0) else 0.0001))
    yp = yv - (((yv - yk) * (zv - zp)) / ((zv - zk) if (zv - zk != 0) else 0.0001))

    return [int(res * (xp * 50 + 100)), int(res * (yp * 50 + 100))]


def calculateSurfaceNormal(triangle):
    U = substract(triangle[1], triangle[0])
    V = substract(triangle[2], triangle[0])

    Nx = (U[1] * V[2]) - (U[2] * V[1])
    Ny = (U[2] * V[0]) - (U[0] * V[2])
    Nz = (U[0] * V[1]) - (U[1] * V[0])

    return [Nx / 5, Ny / 5, Nz / 5]


def distance(a, b):
    """Return the distance between two vectors using the pythagorean theory
    """
    # return math.sqrt(math.pow((b[0] - a[0]), 2) + math.pow((b[1] - a[1]), 2) + math.pow((b[2] - a[2]), 2))
    # # print(a,b)
    # # return math.hypot(a, b)
    return (((b[0] - a[0]) ** 2) + ((b[1] - a[1]) ** 2) + ((b[2] - a[2]) ** 2)) ** 0.5


def distance_squared(a, b):
    """Return the distance squaredbetween two vectors using the pythagorean theory
    """
    return ((b[0] - a[0]) ** 2) + ((b[1] - a[1]) ** 2) + ((b[2] - a[2]) ** 2)


def degreesToRadians(angle):
    """Simple degrees-to-radians converter
    """
    return (angle * 2 * math.pi) / 360


def radiansToDegrees(angle):
    """Simple radians-to-degrees converter
    """
    return (angle * 360) / (math.pi * 2)


def angleBetweenVectors(a, b):
    """Returns the angle between two vectors, in degreesToRadians
    """

    # calculates the cosine
    cos = (a[0] * b[0] + a[1] * b[1] + a[2] * b[2]) / (
        distance([0, 0, 0], a) * distance([0, 0, 0], b)
    )

    return radiansToDegrees(math.acos(cos))


def sortFacesByDistance(arr):
    """Uses Insertion Sort to sort all the faces by the distance from the center to the camera
    """

    for i in tqdm(range(1, len(arr)), ascii=True, desc="Sorting faces"):

        key = arr[i]

        j = i - 1
        while j >= 0 and key[0] < arr[j][0]:
            arr[j + 1] = arr[j]
            j -= 1
        arr[j + 1] = key


def colorLerp(color1, color2, t):
    return (
        int(color1[0] + ((color2[0] - color1[0]) * t)),
        int(color1[1] + ((color2[1] - color1[1]) * t)),
        int(color1[2] + ((color2[2] - color1[2]) * t)),
    )


def vectorLerpRound(v1, v2, t):
    return [int(v1[0] + ((v2[0] - v1[0]) * t)), int(v1[1] + ((v2[1] - v1[1]) * t))]


def vectorLerp(v1, v2, t):
    return [
        v1[0] + ((v2[0] - v1[0]) * t),
        v1[1] + ((v2[1] - v1[1]) * t),
        v1[2] + ((v2[2] - v1[2]) * t),
    ]


def draw_line(start, end, startColor, endColor):
    """Bresenham's Line Algorithm
    Produces a list of tuples from start and end
 
    >>> points1 = get_line((0, 0), (3, 4))
    >>> points2 = get_line((3, 4), (0, 0))
    >>> assert(set(points1) == set(points2))
    >>> print points1
    [(0, 0), (1, 1), (1, 2), (2, 3), (3, 4)]
    >>> print points2
    [(3, 4), (2, 3), (1, 2), (1, 1), (0, 0)]
    """
    # Setup initial conditions
    x1, y1 = start
    x2, y2 = end
    dx = x2 - x1
    dy = y2 - y1

    # Determine how steep the line is
    is_steep = abs(dy) > abs(dx)

    # Rotate line
    if is_steep:
        x1, y1 = y1, x1
        x2, y2 = y2, x2

    # Swap start and end points if necessary and store swap state
    swapped = False
    if x1 > x2:
        x1, x2 = x2, x1
        y1, y2 = y2, y1
        swapped = True

    # Recalculate differentials
    dx = x2 - x1
    dy = y2 - y1

    # Calculate error
    error = int(dx / 2.0)
    ystep = 1 if y1 < y2 else -1

    # Iterate over bounding box generating points between start and end
    y = y1
    points = []
    for x in range(x1, x2 + 1):
        coord = (y, x) if is_steep else (x, y)
        points.append(coord)
        error -= abs(dy)
        if error < 0:
            y += ystep
            error += dx

    # Reverse the list if the coordinates were swapped
    if swapped:
        points.reverse()

    imgpixels = img.load()
    for pixel in points:
        dist = distance([end[0], end[1], 0], [start[0], start[1], 0])
        if dist != 0:
            try:
                imgpixels[pixel[0], pixel[1]] = colorLerp(
                    startColor,
                    endColor,
                    distance([pixel[0], pixel[1], 0], [start[0], start[1], 0]) / dist,
                )
            except IndexError:
                noname = "wut"
        else:
            try:
                imgpixels[pixel[0], pixel[1]] = endColor
            except IndexError:
                noname = "wut"


animationScript = settings.animationScript

shadingMode = settings.shadingMode

yRotation = settings.yRotation
xRotation = settings.xRotation
zRotation = settings.zRotation

res = settings.imageResolution / 100

curves = settings.curves

for crv in range(len(curves)):
    curves[crv][0][1] = 1 - curves[crv][0][1]
    curves[crv][1][1] = 1 - curves[crv][1][1]
    curves[crv][2][1] = 1 - curves[crv][2][1]
    curves[crv][3][1] = 1 - curves[crv][3][1]

diffuseColor = settings.diffuseColor
gouraudColor = settings.gouraudColor
phongColor = settings.gouraudColor


meshPosition = settings.meshPosition

renderFaceNormals = settings.renderFaceNormals
renderHardVertexNormals = settings.renderHardVertexNormals
renderSoftVertexNormals = settings.renderSoftVertexNormals
renderMeshCenter = settings.renderMeshCenter
renderFaceCenters = settings.renderFaceCenters
renderLightSource = settings.renderLightSource

cameraPosition = settings.cameraPosition

imagePlaneZ = settings.imagePlaneZ

animationScript = animationScript.replace(" ", "")
animationScript = animationScript.split("/")
renderNumber = 0
for key in animationScript:
    key = key.split(";")

    repeatMode = "steps"
    interpolation = "bezier"
    bezierNumber = 2
    numOfTimes = 1
    for item in key:
        if item.endswith("f"):
            repeatMode = "steps"
            numOfTimes = int(item[:-1])

    xRotationTemp = xRotation
    yRotationTemp = yRotation
    zRotationTemp = zRotation
    meshPostion0Temp = meshPosition[0]
    meshPostion1Temp = meshPosition[1]
    meshPostion2Temp = meshPosition[2]
    resTemp = res
    for current in range(numOfTimes):
        renderNumber += 1
        time = (current + 1) / numOfTimes
        for item in key:
            item = item.split(",")
            for att in item:
                if att == "t":
                    repeatMode = "times"
                if att == "s":
                    repeatMode = "steps"
                if att == "i=l":
                    interpolation = "linear"
                if att.startswith("i=b"):
                    interpolation = "bezier"
                    bezierNumber = int(att[3:])
            for att in item:
                pointy = (
                    1
                    - intersection.computeIntersections(
                        [
                            curves[bezierNumber][0][0],
                            curves[bezierNumber][1][0],
                            curves[bezierNumber][2][0],
                            curves[bezierNumber][3][0],
                        ],
                        [
                            curves[bezierNumber][0][1],
                            curves[bezierNumber][1][1],
                            curves[bezierNumber][2][1],
                            curves[bezierNumber][3][1],
                        ],
                        [time, time],
                        [-1, 2],
                    )[0][1]
                )
                if "+" in att:
                    att = att.split("+")
                    if att[0] == "xRotation":
                        xRotation = (
                            (
                                (xRotation + float(att[1]))
                                if repeatMode == "times"
                                else (xRotation + float(att[1]) / numOfTimes)
                            )
                            if (interpolation == "linear")
                            else (
                                (xRotationTemp + float(att[1]) * numOfTimes * pointy)
                                if repeatMode == "times"
                                else (xRotationTemp + float(att[1]) * pointy)
                            )
                        )
                    if att[0] == "yRotation":
                        yRotation = (
                            (
                                (yRotation + float(att[1]))
                                if repeatMode == "times"
                                else (yRotation + float(att[1]) / numOfTimes)
                            )
                            if (interpolation == "linear")
                            else (
                                (yRotationTemp + float(att[1]) * numOfTimes * pointy)
                                if repeatMode == "times"
                                else (yRotationTemp + float(att[1]) * pointy)
                            )
                        )
                    if att[0] == "zRotation":
                        zRotation = (
                            (
                                (zRotation + float(att[1]))
                                if repeatMode == "times"
                                else (zRotation + float(att[1]) / numOfTimes)
                            )
                            if (interpolation == "linear")
                            else (
                                (zRotationTemp + float(att[1]) * numOfTimes * pointy)
                                if repeatMode == "times"
                                else (zRotationTemp + float(att[1]) * pointy)
                            )
                        )
                    if att[0] == "meshPositionX":
                        meshPosition[0] = (
                            (
                                (meshPosition[0] + float(att[1]))
                                if repeatMode == "times"
                                else (meshPosition[0] + float(att[1]) / numOfTimes)
                            )
                            if (interpolation == "linear")
                            else (
                                (meshPostion0Temp + float(att[1]) * numOfTimes * pointy)
                                if repeatMode == "times"
                                else (meshPostion0Temp + float(att[1]) * pointy)
                            )
                        )
                    if att[0] == "meshPositionY":
                        meshPosition[1] = (
                            (
                                (meshPosition[1] + float(att[1]))
                                if repeatMode == "times"
                                else (meshPosition[1] + float(att[1]) / numOfTimes)
                            )
                            if (interpolation == "linear")
                            else (
                                (meshPostion1Temp + float(att[1]) * numOfTimes * pointy)
                                if repeatMode == "times"
                                else (meshPostion1Temp + float(att[1]) * pointy)
                            )
                        )
                    if att[0] == "meshPositionZ":
                        meshPosition[2] = (
                            (
                                (meshPosition[2] + float(att[1]))
                                if repeatMode == "times"
                                else (meshPosition[2] + float(att[1]) / numOfTimes)
                            )
                            if (interpolation == "linear")
                            else (
                                (meshPostion2Temp + float(att[1]) * numOfTimes * pointy)
                                if repeatMode == "times"
                                else (meshPostion2Temp + float(att[1]) * pointy)
                            )
                        )
                    if att[0] == "resolution":
                        res = (
                            (
                                (res + float(att[1]) / 100)
                                if repeatMode == "times"
                                else (res + (float(att[1]) / 100) / numOfTimes)
                            )
                            if (interpolation == "linear")
                            else (
                                (resTemp + (float(att[1]) / 100) * numOfTimes * pointy)
                                if repeatMode == "times"
                                else (resTemp + (float(att[1]) / 100) * pointy)
                            )
                        )
                if "-" in att and "=" not in att:
                    att = att.split("-")
                    if att[0] == "xRotation":
                        xRotation = (
                            (
                                (xRotation - float(att[1]))
                                if repeatMode == "times"
                                else (xRotation - float(att[1]) / numOfTimes)
                            )
                            if (interpolation == "linear")
                            else (
                                (xRotationTemp - float(att[1]) * numOfTimes * pointy)
                                if repeatMode == "times"
                                else (xRotationTemp - float(att[1]) * pointy)
                            )
                        )
                    if att[0] == "yRotation":
                        yRotation = (
                            (
                                (yRotation - float(att[1]))
                                if repeatMode == "times"
                                else (yRotation - float(att[1]) / numOfTimes)
                            )
                            if (interpolation == "linear")
                            else (
                                (yRotationTemp - float(att[1]) * numOfTimes * pointy)
                                if repeatMode == "times"
                                else (yRotationTemp - float(att[1]) * pointy)
                            )
                        )
                    if att[0] == "zRotation":
                        zRotation = (
                            (
                                (zRotation - float(att[1]))
                                if repeatMode == "times"
                                else (zRotation - float(att[1]) / numOfTimes)
                            )
                            if (interpolation == "linear")
                            else (
                                (zRotationTemp - float(att[1]) * numOfTimes * pointy)
                                if repeatMode == "times"
                                else (zRotationTemp - float(att[1]) * pointy)
                            )
                        )
                    if att[0] == "meshPositionX":
                        meshPosition[0] = (
                            (
                                (meshPosition[0] - float(att[1]))
                                if repeatMode == "times"
                                else (meshPosition[0] - float(att[1]) / numOfTimes)
                            )
                            if (interpolation == "linear")
                            else (
                                (meshPostion0Temp - float(att[1]) * numOfTimes * pointy)
                                if repeatMode == "times"
                                else (meshPostion0Temp - float(att[1]) * pointy)
                            )
                        )
                    if att[0] == "meshPositionY":
                        meshPosition[1] = (
                            (
                                (meshPosition[1] - float(att[1]))
                                if repeatMode == "times"
                                else (meshPosition[1] - float(att[1]) / numOfTimes)
                            )
                            if (interpolation == "linear")
                            else (
                                (meshPostion1Temp - float(att[1]) * numOfTimes * pointy)
                                if repeatMode == "times"
                                else (meshPostion1Temp - float(att[1]) * pointy)
                            )
                        )
                    if att[0] == "meshPositionZ":
                        meshPosition[2] = (
                            (
                                (meshPosition[2] - float(att[1]))
                                if repeatMode == "times"
                                else (meshPosition[2] - float(att[1]) / numOfTimes)
                            )
                            if (interpolation == "linear")
                            else (
                                (meshPostion2Temp - float(att[1]) * numOfTimes * pointy)
                                if repeatMode == "times"
                                else (meshPostion2Temp - float(att[1]) * pointy)
                            )
                        )
                    if att[0] == "resolution":
                        res = (
                            (
                                (res - float(att[1]) / 100)
                                if repeatMode == "times"
                                else (res - (float(att[1]) / 100) / numOfTimes)
                            )
                            if (interpolation == "linear")
                            else (
                                (resTemp - (float(att[1]) / 100) * numOfTimes * pointy)
                                if repeatMode == "times"
                                else (resTemp - (float(att[1]) / 100) * pointy)
                            )
                        )
                if "=" in att:
                    att = att.split("=")
                    if att[0] == "xRotation":
                        xRotation = (
                            ((xRotation + (float(att[1]) - xRotationTemp)) / numOfTimes)
                            if (interpolation == "linear")
                            else (
                                xRotationTemp + (float(att[1]) - xRotationTemp) * pointy
                            )
                        )
                    if att[0] == "yRotation":
                        yRotation = (
                            ((yRotation + (float(att[1]) - yRotationTemp)) / numOfTimes)
                            if (interpolation == "linear")
                            else (
                                yRotationTemp + (float(att[1]) - yRotationTemp) * pointy
                            )
                        )
                    if att[0] == "zRotation":
                        zRotation = (
                            ((zRotation + (float(att[1]) - zRotationTemp)) / numOfTimes)
                            if (interpolation == "linear")
                            else (
                                zRotationTemp + (float(att[1]) - zRotationTemp) * pointy
                            )
                        )
                    if att[0] == "meshPositionX":
                        meshPosition[0] = (
                            (
                                (meshPosition[0] + (float(att[1]) - meshPostion0Temp))
                                / numOfTimes
                            )
                            if (interpolation == "linear")
                            else (
                                meshPostion0Temp
                                + (float(att[1]) - meshPostion0Temp) * pointy
                            )
                        )
                    if att[0] == "meshPositionY":
                        meshPosition[1] = (
                            (
                                (meshPosition[1] + (float(att[1]) - meshPostion1Temp))
                                / numOfTimes
                            )
                            if (interpolation == "linear")
                            else (
                                meshPostion1Temp
                                + (float(att[1]) - meshPostion1Temp) * pointy
                            )
                        )
                    if att[0] == "meshPositionZ":
                        meshPosition[2] = (
                            (
                                (meshPosition[2] + (float(att[1]) - meshPostion2Temp))
                                / numOfTimes
                            )
                            if (interpolation == "linear")
                            else (
                                meshPostion2Temp
                                + (float(att[1]) - meshPostion2Temp) * pointy
                            )
                        )
                    if att[0] == "resolution":
                        res = (
                            ((res + (float(att[1]) / 100 - resTemp)) / numOfTimes)
                            if (interpolation == "linear")
                            else (resTemp + (float(att[1]) / 100 - resTemp) * pointy)
                        )

        initialVertices = []

        verticesOnScreen = []
        rotatedVerticesXYZ = []
        for vertex in tqdm(meshes.meshVertices, ascii=True, desc="Loading vertices"):
            initialVertices.append(vertex)

        verticesPosition = list(map(lambda x: add(x, meshPosition), initialVertices))
        centerPoint = center(verticesPosition)

        lightVector = settings.lightVector
        lightDirectionX = settings.lightDirectionX
        lightDirectionY = settings.lightDirectionY
        lightDirectionZ = settings.lightDirectionZ

        lightDirectionX = degreesToRadians(lightDirectionX)
        lightDirectionY = degreesToRadians(lightDirectionY)
        lightDirectionZ = degreesToRadians(lightDirectionZ)

        hardVertexNormals = []
        softVertexNormals = []

        rotatedVerticesXYZ = rotatePoints(
            verticesPosition,
            centerPoint,
            degreesToRadians(yRotation),
            degreesToRadians(xRotation),
            degreesToRadians(zRotation),
        )

        lightVector = rotatePoints(
            lightVector, [0, 0, 0], lightDirectionY, lightDirectionX, lightDirectionZ
        )

        newLightPoint = calculateScreenPoint(lightVector[0])

        for v in tqdm(rotatedVerticesXYZ, ascii=True, desc="Calculating screen points"):

            newPoint = calculateScreenPoint(v)

            verticesOnScreen.append(newPoint)

        newCenterPoint = calculateScreenPoint(centerPoint)

        meshFaces = []
        mesh3DFaces = []

        for face in tqdm(meshes.meshFaces, ascii=True, desc="Loading faces"):
            meshFaces.append(
                [
                    verticesOnScreen[face[0]],
                    verticesOnScreen[face[1]],
                    verticesOnScreen[face[2]],
                ]
            )
            mesh3DFaces.append(
                [
                    rotatedVerticesXYZ[face[0]],
                    rotatedVerticesXYZ[face[1]],
                    rotatedVerticesXYZ[face[2]],
                ]
            )

        gridCenterOnScreen = calculateScreenPoint([0, 0, 0])

        faceNormalsOnScreen = []
        faceCenters = []

        for face in tqdm(mesh3DFaces, ascii=True, desc="Calculating face centers"):
            faceCenters.append(center([face[0], face[1], face[2]]))

        sortedDistances = [
            [distance_squared(cen, cameraPosition), c]
            for c, cen in enumerate(faceCenters)
        ]
        # print(sortedDistances)
        # sortFacesByDistance(sortedDistances)
        # # TODO if you reverse sortFacesByDistance anyways, why not sort them in reverse?
        # sortedDistances.reverse()
        sortedDistances = sorted(sortedDistances, key=lambda d: d[0], reverse=True)
        tempMeshFaces = []
        tempMesh3DFaces = []
        for d in sortedDistances:
            tempMeshFaces.append(meshFaces[d[1]])
            tempMesh3DFaces.append(mesh3DFaces[d[1]])
        meshFaces = tempMeshFaces
        mesh3DFaces = tempMesh3DFaces

        k = 0
        faceNormals = []

        for face in tqdm(mesh3DFaces, ascii=True, desc="Calculating face normals"):
            faceNormal = calculateSurfaceNormal(face)

            faceNormals.append(
                [
                    faceNormal[0] + faceCenters[k][0],
                    faceNormal[1] + faceCenters[k][1],
                    faceNormal[2] + faceCenters[k][2],
                ]
            )

            k += 1

        for vertex in rotatedVerticesXYZ:
            hardVertexNormals.append([])

        if settings.calculateVertexNormals:
            h = 0
            for face in tqdm(
                mesh3DFaces, ascii=True, desc="Calculating vertex normals"
            ):
                vertexNormals = []

                for vertex in face:
                    hardVertexNormals[rotatedVerticesXYZ.index(vertex)].append(
                        [
                            faceNormals[h][0] - faceCenters[h][0] + vertex[0],
                            faceNormals[h][1] - faceCenters[h][1] + vertex[1],
                            faceNormals[h][2] - faceCenters[h][2] + vertex[2],
                        ]
                    )
                h += 1

        for normals in hardVertexNormals:
            softNormal = center(normals)
            softVertexNormals.append(softNormal)

        if renderFaceNormals:
            for n in faceNormals:
                newPoint = calculateScreenPoint(n)
                faceNormalsOnScreen.append(newPoint)

        faceCentersScreen = []
        if renderFaceCenters:
            for n in faceCenters:
                newPoint = calculateScreenPoint(n)
                faceCentersScreen.append(newPoint)

        hardVerticesNormalsOnScreen = []
        softVerticesNormalsOnScreen = []

        if renderHardVertexNormals:
            for v in hardVertexNormals:
                vertexNormalsOnScreen = []
                for n in v:
                    newPoint = calculateScreenPoint(n)
                    vertexNormalsOnScreen.append(newPoint)
                hardVerticesNormalsOnScreen.append(vertexNormalsOnScreen)

        if renderSoftVertexNormals:
            for n in softVertexNormals:
                newPoint = calculateScreenPoint(n)
                softVerticesNormalsOnScreen.append(newPoint)

        facesColors = []
        vertexColors = []

        j = 0
        for normal in tqdm(faceNormals, ascii=True, desc="Calculating face colors"):
            lightNormalAngle = angleBetweenVectors(
                [
                    normal[0] - faceCenters[j][0],
                    normal[1] - faceCenters[j][1],
                    normal[2] - faceCenters[j][2],
                ],
                lightVector[0],
            )

            if lightNormalAngle <= 90:
                if (
                    shadingMode == settings.ShadingMode.flatDiffuse
                    or shadingMode == settings.ShadingMode.gouraud
                ):
                    facesColors.append((2, 2, 2))
                elif shadingMode == settings.ShadingMode.unlit:
                    facesColors.append(
                        (
                            settings.unlitColor[0],
                            settings.unlitColor[1],
                            settings.unlitColor[2],
                        )
                    )
            else:
                if (
                    shadingMode == settings.ShadingMode.flatDiffuse
                    or shadingMode == settings.ShadingMode.gouraud
                ):
                    lightNormalAngle -= 90
                    color = (
                        int(lightNormalAngle * diffuseColor[0] / 90),
                        int(lightNormalAngle * diffuseColor[1] / 90),
                        int(lightNormalAngle * diffuseColor[2] / 90),
                    )
                    facesColors.append(color)
                elif shadingMode == settings.ShadingMode.unlit:
                    facesColors.append(
                        (
                            settings.unlitColor[0],
                            settings.unlitColor[1],
                            settings.unlitColor[2],
                        )
                    )
            j += 1

        j = 0
        for normal in tqdm(
            softVertexNormals, ascii=True, desc="Calculating vertex colors"
        ):
            lightNormalAngle = angleBetweenVectors(
                [
                    normal[0] - rotatedVerticesXYZ[j][0],
                    normal[1] - rotatedVerticesXYZ[j][1],
                    normal[2] - rotatedVerticesXYZ[j][2],
                ],
                lightVector[0],
            )

            if lightNormalAngle <= 90:
                if shadingMode == settings.ShadingMode.gouraud:
                    vertexColors.append((2, 2, 2))
            else:
                if shadingMode == settings.ShadingMode.gouraud:
                    lightNormalAngle -= 90
                    color = (
                        int(lightNormalAngle * gouraudColor[0] / 90),
                        int(lightNormalAngle * gouraudColor[1] / 90),
                        int(lightNormalAngle * gouraudColor[2] / 90),
                    )
                    vertexColors.append(color)
            j += 1

        width = int(300 * res)
        height = int(300 * res)
        img = Image.new("RGB", (width, height), color="black")

        draw = ImageDraw.Draw(img)

        facesDrawn = []

        s = 0
        for face in mesh3DFaces:
            if (
                angleBetweenVectors(
                    [
                        faceNormals[s][0] - faceCenters[s][0],
                        faceNormals[s][1] - faceCenters[s][1],
                        faceNormals[s][2] - faceCenters[s][2],
                    ],
                    [
                        cameraPosition[0] - faceCenters[s][0],
                        cameraPosition[1] - faceCenters[s][1],
                        cameraPosition[2] - faceCenters[s][2],
                    ],
                )
                < 90
            ):
                facesDrawn.append(1)
            else:
                facesDrawn.append(0)
            s += 1

        i = 0
        for face in tqdm(meshFaces, ascii=True, desc="Rendering faces"):
            if (
                facesDrawn[i] == 1 and settings.backfaceCulling == True
            ) or settings.backfaceCulling == False:
                if (
                    shadingMode == settings.ShadingMode.unlit
                    or shadingMode == settings.ShadingMode.flatDiffuse
                ):
                    draw.polygon(
                        [coord for vertex in face for coord in vertex],
                        fill=facesColors[i],
                    )
                elif shadingMode == settings.ShadingMode.wireframe:
                    draw.polygon(
                        [coord for vertex in face for coord in vertex],
                        outline=(
                            settings.wireframeColor[0],
                            settings.wireframeColor[1],
                            settings.wireframeColor[2],
                        ),
                        fill=(0, 0, 0),
                    )
                elif shadingMode == settings.ShadingMode.gouraud:
                    b = 0
                    dist1 = distance_squared(
                        [face[1][0], face[1][1], 0], [face[0][0], face[0][1], 0]
                    )
                    dist2 = distance_squared(
                        [face[1][0], face[1][1], 0], [face[2][0], face[2][1], 0]
                    )
                    dist3 = distance_squared(
                        [face[0][0], face[0][1], 0], [face[2][0], face[2][1], 0]
                    )
                    if dist1 > dist2:
                        while b < 1:
                            draw_line(
                                vectorLerpRound(face[1], face[0], b),
                                vectorLerpRound(face[1], face[2], b),
                                colorLerp(
                                    vertexColors[
                                        rotatedVerticesXYZ.index(mesh3DFaces[i][1])
                                    ],
                                    vertexColors[
                                        rotatedVerticesXYZ.index(mesh3DFaces[i][0])
                                    ],
                                    b,
                                ),
                                colorLerp(
                                    vertexColors[
                                        rotatedVerticesXYZ.index(mesh3DFaces[i][1])
                                    ],
                                    vertexColors[
                                        rotatedVerticesXYZ.index(mesh3DFaces[i][2])
                                    ],
                                    b,
                                ),
                            )
                            try:
                                b += 1 / int(dist1 ** 0.5)
                            except ZeroDivisionError:
                                b += 1
                    else:
                        while b < 1:
                            draw_line(
                                vectorLerpRound(face[1], face[0], b),
                                vectorLerpRound(face[1], face[2], b),
                                colorLerp(
                                    vertexColors[
                                        rotatedVerticesXYZ.index(mesh3DFaces[i][1])
                                    ],
                                    vertexColors[
                                        rotatedVerticesXYZ.index(mesh3DFaces[i][0])
                                    ],
                                    b,
                                ),
                                colorLerp(
                                    vertexColors[
                                        rotatedVerticesXYZ.index(mesh3DFaces[i][1])
                                    ],
                                    vertexColors[
                                        rotatedVerticesXYZ.index(mesh3DFaces[i][2])
                                    ],
                                    b,
                                ),
                            )
                            try:
                                b += (1 / int(dist2 ** 0.5)) / 2
                            except ZeroDivisionError:
                                b += 1

                    b = 0
                    if dist1 > dist3:
                        while b < 1:
                            draw_line(
                                vectorLerpRound(face[0], face[1], b),
                                vectorLerpRound(face[0], face[2], b),
                                colorLerp(
                                    vertexColors[
                                        rotatedVerticesXYZ.index(mesh3DFaces[i][0])
                                    ],
                                    vertexColors[
                                        rotatedVerticesXYZ.index(mesh3DFaces[i][1])
                                    ],
                                    b,
                                ),
                                colorLerp(
                                    vertexColors[
                                        rotatedVerticesXYZ.index(mesh3DFaces[i][0])
                                    ],
                                    vertexColors[
                                        rotatedVerticesXYZ.index(mesh3DFaces[i][2])
                                    ],
                                    b,
                                ),
                            )
                            try:
                                b += 1 / int(dist1 ** 0.5)
                            except ZeroDivisionError:
                                b += 1
                    else:
                        while b < 1:
                            draw_line(
                                vectorLerpRound(face[0], face[1], b),
                                vectorLerpRound(face[0], face[2], b),
                                colorLerp(
                                    vertexColors[
                                        rotatedVerticesXYZ.index(mesh3DFaces[i][0])
                                    ],
                                    vertexColors[
                                        rotatedVerticesXYZ.index(mesh3DFaces[i][1])
                                    ],
                                    b,
                                ),
                                colorLerp(
                                    vertexColors[
                                        rotatedVerticesXYZ.index(mesh3DFaces[i][0])
                                    ],
                                    vertexColors[
                                        rotatedVerticesXYZ.index(mesh3DFaces[i][2])
                                    ],
                                    b,
                                ),
                            )
                            try:
                                b += (1 / int(dist3 ** 0.5)) / 2
                            except ZeroDivisionError:
                                b += 1
                    b = 0
                    if dist2 > dist3:
                        while b < 1:
                            draw_line(
                                vectorLerpRound(face[2], face[1], b),
                                vectorLerpRound(face[2], face[0], b),
                                colorLerp(
                                    vertexColors[
                                        rotatedVerticesXYZ.index(mesh3DFaces[i][2])
                                    ],
                                    vertexColors[
                                        rotatedVerticesXYZ.index(mesh3DFaces[i][1])
                                    ],
                                    b,
                                ),
                                colorLerp(
                                    vertexColors[
                                        rotatedVerticesXYZ.index(mesh3DFaces[i][2])
                                    ],
                                    vertexColors[
                                        rotatedVerticesXYZ.index(mesh3DFaces[i][0])
                                    ],
                                    b,
                                ),
                            )
                            try:
                                b += 1 / int(dist2 ** 0.5)
                            except ZeroDivisionError:
                                b += 1
                    else:
                        while b < 1:
                            draw_line(
                                vectorLerpRound(face[2], face[1], b),
                                vectorLerpRound(face[2], face[0], b),
                                colorLerp(
                                    vertexColors[
                                        rotatedVerticesXYZ.index(mesh3DFaces[i][2])
                                    ],
                                    vertexColors[
                                        rotatedVerticesXYZ.index(mesh3DFaces[i][1])
                                    ],
                                    b,
                                ),
                                colorLerp(
                                    vertexColors[
                                        rotatedVerticesXYZ.index(mesh3DFaces[i][2])
                                    ],
                                    vertexColors[
                                        rotatedVerticesXYZ.index(mesh3DFaces[i][0])
                                    ],
                                    b,
                                ),
                            )
                            try:
                                b += (1 / int(dist3 ** 0.5)) / 2
                            except ZeroDivisionError:
                                b += 1
                elif shadingMode == settings.ShadingMode.phong:
                    b = 0
                    o = 0
                    normals10 = []
                    normals12 = []
                    normals02 = []

                    vertices10 = []
                    vertices12 = []
                    vertices02 = []

                    colors10 = []
                    colors12 = []
                    colors02 = []
                    dist1 = distance_squared(
                        [face[1][0], face[1][1], 0], [face[0][0], face[0][1], 0]
                    )
                    dist2 = distance_squared(
                        [face[1][0], face[1][1], 0], [face[2][0], face[2][1], 0]
                    )
                    dist3 = distance_squared(
                        [face[0][0], face[0][1], 0], [face[2][0], face[2][1], 0]
                    )

                    while o < 1:
                        normals10.append(
                            vectorLerp(
                                softVertexNormals[
                                    rotatedVerticesXYZ.index(mesh3DFaces[i][1])
                                ],
                                softVertexNormals[
                                    rotatedVerticesXYZ.index(mesh3DFaces[i][0])
                                ],
                                o,
                            )
                        )
                        vertices10.append(
                            vectorLerp(
                                rotatedVerticesXYZ[
                                    rotatedVerticesXYZ.index(mesh3DFaces[i][1])
                                ],
                                rotatedVerticesXYZ[
                                    rotatedVerticesXYZ.index(mesh3DFaces[i][0])
                                ],
                                o,
                            )
                        )

                        try:
                            o += 1 / int(dist1 ** 0.5)
                        except ZeroDivisionError:
                            o += 1
                    o = 0

                    while o < 1:
                        normals12.append(
                            vectorLerp(
                                softVertexNormals[
                                    rotatedVerticesXYZ.index(mesh3DFaces[i][1])
                                ],
                                softVertexNormals[
                                    rotatedVerticesXYZ.index(mesh3DFaces[i][2])
                                ],
                                o,
                            )
                        )
                        vertices12.append(
                            vectorLerp(
                                rotatedVerticesXYZ[
                                    rotatedVerticesXYZ.index(mesh3DFaces[i][1])
                                ],
                                rotatedVerticesXYZ[
                                    rotatedVerticesXYZ.index(mesh3DFaces[i][2])
                                ],
                                o,
                            )
                        )

                        try:
                            o += 1 / int(dist2 ** 0.5)
                        except ZeroDivisionError:
                            o += 1
                    o = 0

                    while o < 1:
                        normals02.append(
                            vectorLerp(
                                softVertexNormals[
                                    rotatedVerticesXYZ.index(mesh3DFaces[i][0])
                                ],
                                softVertexNormals[
                                    rotatedVerticesXYZ.index(mesh3DFaces[i][2])
                                ],
                                o,
                            )
                        )
                        vertices02.append(
                            vectorLerp(
                                rotatedVerticesXYZ[
                                    rotatedVerticesXYZ.index(mesh3DFaces[i][0])
                                ],
                                rotatedVerticesXYZ[
                                    rotatedVerticesXYZ.index(mesh3DFaces[i][2])
                                ],
                                o,
                            )
                        )

                        try:
                            o += 1 / int(dist3 ** 0.5)
                        except ZeroDivisionError:
                            o += 1
                    o = 0
                    for norm in normals10:
                        lightNormalAngle = angleBetweenVectors(
                            [
                                norm[0] - vertices10[normals10.index(norm)][0],
                                norm[1] - vertices10[normals10.index(norm)][1],
                                norm[2] - vertices10[normals10.index(norm)][2],
                            ],
                            lightVector[0],
                        )

                        if lightNormalAngle <= 90:
                            colors10.append((2, 2, 2))
                        else:
                            lightNormalAngle -= 90
                            color = (
                                int(lightNormalAngle * phongColor[0] / 90),
                                int(lightNormalAngle * phongColor[1] / 90),
                                int(lightNormalAngle * phongColor[2] / 90),
                            )
                            colors10.append(color)

                    for norm in normals12:
                        lightNormalAngle = angleBetweenVectors(
                            [
                                norm[0] - vertices12[normals12.index(norm)][0],
                                norm[1] - vertices12[normals12.index(norm)][1],
                                norm[2] - vertices12[normals12.index(norm)][2],
                            ],
                            lightVector[0],
                        )

                        if lightNormalAngle <= 90:
                            colors12.append((2, 2, 2))
                        else:
                            lightNormalAngle -= 90
                            color = (
                                int(lightNormalAngle * phongColor[0] / 90),
                                int(lightNormalAngle * phongColor[1] / 90),
                                int(lightNormalAngle * phongColor[2] / 90),
                            )
                            colors12.append(color)

                    for norm in normals02:
                        lightNormalAngle = angleBetweenVectors(
                            [
                                norm[0] - vertices02[normals02.index(norm)][0],
                                norm[1] - vertices02[normals02.index(norm)][1],
                                norm[2] - vertices02[normals02.index(norm)][2],
                            ],
                            lightVector[0],
                        )

                        if lightNormalAngle <= 90:
                            colors02.append((2, 2, 2))
                        else:
                            lightNormalAngle -= 90
                            color = (
                                int(lightNormalAngle * phongColor[0] / 90),
                                int(lightNormalAngle * phongColor[1] / 90),
                                int(lightNormalAngle * phongColor[2] / 90),
                            )
                            colors02.append(color)

                    try:
                        if dist1 > dist2:
                            while b < 1:
                                draw_line(
                                    vectorLerpRound(face[1], face[0], b),
                                    vectorLerpRound(face[1], face[2], b),
                                    colors10[int(b * int(dist1 ** 0.5)) - 1],
                                    colors12[int(b * int(dist2 ** 0.5)) - 1],
                                )
                                try:
                                    b += 1 / int(dist1 ** 0.5)
                                except ZeroDivisionError:
                                    b += 1
                        else:
                            while b < 1:
                                draw_line(
                                    vectorLerpRound(face[1], face[0], b),
                                    vectorLerpRound(face[1], face[2], b),
                                    colors10[int(b * int(dist1 ** 0.5)) - 1],
                                    colors12[int(b * int(dist2 ** 0.5)) - 1],
                                )
                                try:
                                    b += (1 / int(dist2 ** 0.5)) / 2
                                except ZeroDivisionError:
                                    b += 1
                    except IndexError:
                        pass

                    b = 0
                    try:
                        if dist1 > dist3:
                            while b < 1:
                                draw_line(
                                    vectorLerpRound(face[0], face[1], b),
                                    vectorLerpRound(face[0], face[2], b),
                                    colors10[int(1 - b * int(dist1 ** 0.5)) - 1],
                                    colors02[int(b * int(dist3 ** 0.5)) - 1],
                                )
                                try:
                                    b += 1 / int(dist1 ** 0.5)
                                except ZeroDivisionError:
                                    b += 1
                        else:
                            while b < 1:
                                draw_line(
                                    vectorLerpRound(face[0], face[1], b),
                                    vectorLerpRound(face[0], face[2], b),
                                    colors10[int(1 - b * int(dist1 ** 0.5)) - 1],
                                    colors02[int(b * int(dist3 ** 0.5)) - 1],
                                )
                                try:
                                    b += (1 / int(dist3 ** 0.5)) / 2
                                except ZeroDivisionError:
                                    b += 1
                    except IndexError:
                        pass
                    b = 0
                    try:
                        if dist2 > dist3:
                            while b < 1:
                                draw_line(
                                    vectorLerpRound(face[2], face[1], b),
                                    vectorLerpRound(face[2], face[0], b),
                                    colors12[int(1 - b * int(dist2 ** 0.5)) - 1],
                                    colors02[int(1 - b * int(dist3 ** 0.5)) - 1],
                                )
                                try:
                                    b += 1 / int(dist2 ** 0.5)
                                except ZeroDivisionError:
                                    b += 1
                        else:
                            while b < 1:
                                draw_line(
                                    vectorLerpRound(face[2], face[1], b),
                                    vectorLerpRound(face[2], face[0], b),
                                    colors12[int(1 - b * int(dist2 ** 0.5)) - 1],
                                    colors02[int(1 - b * int(dist3 ** 0.5)) - 1],
                                )
                                try:
                                    b += (1 / int(dist3 ** 0.5)) / 2
                                except ZeroDivisionError:
                                    b += 1
                    except IndexError:
                        pass
            if renderFaceNormals:
                draw.line(
                    [
                        faceCentersScreen[i][0],
                        faceCentersScreen[i][1],
                        faceNormalsOnScreen[i][0],
                        faceNormalsOnScreen[i][1],
                    ],
                    (160, 0, 0),
                )
            i += 1

        if renderHardVertexNormals:
            m = 0
            for vertex in hardVerticesNormalsOnScreen:
                for normal in vertex:
                    draw.line(
                        [
                            verticesOnScreen[m][0],
                            verticesOnScreen[m][1],
                            normal[0],
                            normal[1],
                        ],
                        (0, 0, 200),
                    )
                m += 1

        if renderSoftVertexNormals:
            z = 0
            for normal in softVerticesNormalsOnScreen:
                try:
                    draw.line(
                        [
                            verticesOnScreen[z][0],
                            verticesOnScreen[z][1],
                            normal[0],
                            normal[1],
                        ],
                        vertexColors[z],
                    )
                except IndexError:
                    noname = "wut"
                z += 1

        if renderFaceCenters:
            for point in faceCentersScreen:
                draw.ellipse(
                    (point[0] - 1, point[1] - 1, point[0] + 1, point[1] + 1), fill="red"
                )

        if renderMeshCenter:
            draw.ellipse(
                (
                    newCenterPoint[0] - 1,
                    newCenterPoint[1] - 1,
                    newCenterPoint[0] + 1,
                    newCenterPoint[1] + 1,
                ),
                fill="gray",
            )

        if renderLightSource:
            draw.ellipse(
                (
                    gridCenterOnScreen[0] - 2,
                    gridCenterOnScreen[1] - 2,
                    gridCenterOnScreen[0] + 2,
                    gridCenterOnScreen[1] + 2,
                ),
                fill="lightgray",
            )
            draw.line(
                [
                    gridCenterOnScreen[0],
                    gridCenterOnScreen[1],
                    newLightPoint[0],
                    newLightPoint[1],
                ],
                "yellow",
            )

        img.save("render" + str(renderNumber) + ".png")
img.show()
