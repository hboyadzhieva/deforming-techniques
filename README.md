# deforming-techniques
Implementation of 2D Grid Deformation and Free-Form Deformation techniques with python

## About
Algorithms used for warping/deforming objects.

1. **2D grid deformation**  
Create a lattice around a 2D model. Displace the lattice vertices and calculate the new position of each vertex from the model using bilinear interpolation.

2. **Free-form deformation**  
Create a lattice around a 3D model. Displace the lattice vertices and calculate the new position of each vertex from the model using Bezier volumes.

## Run locally
- Clone the project  
`git clone https://github.com/hboyadzhieva/deforming-techniques.git`
- go to project folder and install external dependencies  
`pip install -r requirements.txt`
- run  
`python main.py`
