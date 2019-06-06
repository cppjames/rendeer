filename = input("Please enter the name of the .obj file: ")
objfile = open(filename, "r")

objContent = objfile.read()
objLines = objContent.split("\n")

finalString = "meshVertices = ["
vertices = []
coords = ""

for line in objLines:
    if line[:2] == "v ":
        line = line[2:]
        coords = "[" + line.replace(" ", ", ") + "]"

        finalString += coords + ","
finalString = finalString[:-1]
finalString += "]\n\nmeshFaces = ["

for line in objLines:
    if line[:2] == "f ":
        vertices = []
        line = line[2:]
        elements = line.split(" ")
        for element in elements:
            element = element.split("/")
            vertices.append(int(element[0]) - 1)
        finalString += str(vertices) + ","
finalString = finalString[:-1]
finalString += "]"

meshesFile = open("meshes.py", "w")
meshesFile.write(finalString)
