# Critical
!!NEEDS NUMPY AND NUMBA TO FUNCTION!!

# General Info
Installation:
```python
pip install OskarFieldModule

```
Importing the class:
```
from Field import Field
```

Class for calculating and plotting electromagnetic fields and potentials in 3d. Requires numba and numpy

# Functions

## Efield, Epot
Field and potential from point charge\n
efieldLine, epotLine - Field and potential from line charge parallel to x, y or z axis. Can be placed anywhere in 3D\n


## EfieldCircle, EpotCircle 
Field and potential from circle charge in origin\n

## BfieldLine
Field from line current parallel to x, y or z axis. Can be placed anywhere in 3D\n

## BfieldCircle
Field from circular current in origin\n\n

## PlotVector
Plots vector field. Customize colorscheme, density, and figsize\n

## PlotContour
Plots vector field. Customize colorscheme, levels, norm and figsize\n

## PlotCircle
Plots circle just using radius
    
# See example code below:


## Point Charge Example
```python

L = 2
N = 10
Q = [1.0]
r_Q = np.array([[0.0, 0.0, 0.0]])
plane = 'xy'

rx, ry, Ex, Ey = Field.CalculateEfield(L, N, Q, r_Q, plane)
Field.PlotVector(rx, ry, Ex, Ey, 'quiver', show = True, save=True, 
                 name = os.path.join(os.path.dirname(__file__), "Example_Figures", "PointCharge.pdf"))
N = 100
rx, ry, V = Field.CalculateEpot(L, N, Q, r_Q, plane)
Field.PlotContour(rx, ry, V, show=True, save = True, 
                  name = os.path.join(os.path.dirname(__file__), "Example_Figures", "PointChargeContour.pdf"))
```
<img src="https://github.com/Oskar-Idland/Oskar-Field-Module/blob/main/ELMAG_Module/src/Example_Figures/PointCharge.pdf" alt="Point Charge" title="Point Charge">


![Point Charge](https://github.com/Oskar-Idland/Oskar-Field-Module/blob/main/ELMAG_Module/src/Example_Figures/PointCharge.pdf?raw=true "Point Charge")
## Double Point Charge Example
```python     

L = 2
N = 100
Q = [2.0, -2.0]
r_Q = [-1, 1]
r_Q = np.array([[-1.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
plane = 'xy'

rx, ry, Ex, Ey = Field.CalculateEfield(L, N, Q, r_Q, plane)
Field.PlotVector(rx, ry, Ex, Ey, 'stream', show = True, broken_streamlines = False, save = True, 
                 name = os.path.join(os.path.dirname(__file__), "Example_Figures", "DoublePointCharge.pdf"))
```
## Double Line Charge Example
```python       

L = 2
N = 100
line_charges = [-1, 1]
line_lengths = [1, 1]
line_center_coords = [[0, 0, -1], [0, 0, 1]]
axis = ['x', 'x']
plane = 'xz'

rx, rz, Ex, Ez = Field.CalculateEfieldLine(L, N, line_charges, line_lengths, line_center_coords, axis, plane)
rx, rz, V = Field.CalculateEpotLine(L, N, line_charges, line_lengths, line_center_coords, axis, plane)

Field.PlotVector(rx, rz, Ex, Ez, 'stream', broken_streamlines = False, show = True, save = True, 
                 name = os.path.join(os.path.dirname(__file__), "Example_Figures", "DoubleLineCharge.pdf"))

Field.PlotContour(rx, rz, V, show = True, norm = 'linear', save = True, 
                  name = os.path.join(os.path.dirname(__file__), "Example_Figures", "DoubleLineChargeContour.pdf"))
```
## Circular Charge Example
```python

L = 5
N = 100
circle_charge = [5]
radius = [2]
plane = 'yz'
plane_circles = ['yz']
ry, rz, Ey, Ez = Field.CalculateEfieldCircle(L, N, circle_charge, radius, plane, plane_circles)
Field.PlotVector(ry, rz, Ey, Ez, 'stream', show = False, equal = True)

t = np.linspace(0, 2*np.pi, 100)
plt.plot(radius[0]*np.cos(t), radius[0]*np.sin(t))
plt.show()

N = 500
ry, rz, V = Field.CalculateEpotCircle(L, N, circle_charge, radius, plane, plane_circles)
Field.PlotContour(ry, rz, V, show = False, equal = True)
Field.PlotCircle(radius[0], show = True, save = True, 
                 name = os.path.join(os.path.dirname(__file__), "Example_Figures", "CircularCharge.pdf"))
```
## Line Current Example
```python

L = 5
N = 24
line_currents = [5]
line_lengths = [1]
line_center_coords = [[0.0, 0.0, 0.0]]
axis = ['x']
plane = 'yz'

rx, rz, Bx, Bz = Field.CalculateBfieldLine(L, N, line_currents, line_lengths, line_center_coords, axis, plane)

Field.PlotVector(rx, rz, Bx, Bz, 'quiver', title = 'Magnetic Field from Line Current', show = True, save = True, 
                 name = os.path.join(os.path.dirname(__file__), "Example_Figures", "LineCurrent.pdf"))
```
## Circular Current Example
```python

L = 8
N = 40
circle_currents = [5]
radii = [5]
plane = 'xz'
circle_planes = ['xy']
rx, rz, Bx, Bz = Field.CalculateBfieldCircle(L, N, circle_currents, radii, plane, circle_planes)

Field.PlotVector(rx, rz, Bx, Bz, 'stream', title = 'Magnetic Field from Circular Line Current', broken_streamlines=False, 
                 show = True, cmap = 'inferno', density = .5, save = True, 
                 name = os.path.join(os.path.dirname(__file__), "Example_Figures", "CircularCurrent.pdf"))
```

[Github](https://github.com/Oskar-Idland/Oskar-Field-Module)

