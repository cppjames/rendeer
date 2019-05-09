![alt text](https://i.imgur.com/WJgYoz2.png "Rendeer")

# Rendeer
> A basic **3D renderer** made in Python, using the library Pillow.

![alt text](https://i.imgur.com/PFlAlpJ.png "Two awesome renders")
### Note: This application is currently in development. It has bugs and imperfections.

## Getting started
These instructions will guide you through the installation and usage of Rendeer.

### Requirements
This project was made using Python 3.6. You need to go to [the official Python website](https://www.python.org/downloads/release/python-368/) and install Python 3.6 or a newer version.
You also need the Pillow library, which you can install with the command:
```
pip install pillow
```

### Usage
You should be able to run the program by simply opening `rendeer.py`.
It will take **up to half a minute** for the mesh to be rendered.
The image you will see should be similar (ideally, identical to) [this](https://i.imgur.com/JmAX3Fq.png).

If you want to import a mesh from an .obj file, you now can! Simply open `wavefront_importer.py` and specify the path to the .obj file. 
#### Note: *Please* triangulate your mesh before importing it unless you want it to look like [this](https://i.imgur.com/vbUDdbK.png)
If, for example, you want to import a mesh called `tree.obj` which is in the same folder as `wavefront_importer.py`, you would need to type 
```
tree.obj
```
If it's in a folder called "Trees" that's in the same folder as `wavefront_importer.py`, you would need to type
```
Trees\tree.obj
```
#### Note: The three files (`rendeer.py`, `meshes.py` and `wavefront_importer.py`) have to be all in the same directory for everything to work properly.
The data from the .obj file will be saved in a file called `meshes.py`. You can now run the renderer and it will display your mesh in an image. **You might have to move/rotate your mesh in order to display it as you like.**

## Changelog
### What's new in version 2.3.0?
- shading modes! choose between `flat diffuse`, `unlit` and `wireframe`!
- progress bars! stop staring into the void for a minute straight, now you can watch the program doing things
- variables that the user should modify are in a separate file now (`settings.py`)
- fixed problem with the mesh being too bright where light doesn't touch it

#### v2.2.0
- now it has a Wavefront importer for all your amazing 3D models.
- it sorts faces by distance so it can now render them in the (almost) correct order.
- it can now calculate hard vertex normals (don't know why I made this).
- I've fixed some of the BAD CODE.
