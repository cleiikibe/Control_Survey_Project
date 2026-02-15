# COMPUTATION SYSTEM
This project deals with computation of provisional coordinates and final coordinates of unknown points in triangulation. Provision of horizontal controls by Triangulation Method.
A professional Python implementation for precise coordinate determination in land surveying using classical *intersection/resection* and *cut computation* methods. Designed for surveyors, geomatics engineers, and students requiring rigorous field computation workflows.

## ✨ Key Features

- *Dual computation modules*:
  - *Provisional Coordinates*: Determine unknown station coordinates from two known control points using either direct bearings or observed angles (α/β)
  - *Cut Computations*: Refine coordinates through trigonometric adjustment using cotangent, cosecant, and secant identities
- *Survey-grade precision*:
  - DMS (Degrees-Minutes-Seconds) ↔ Decimal degree conversion with proper rounding
  - 6-decimal precision handling before trigonometric operations (per surveying standards)
  - Robust error checking for collinear points and invalid geometries
- *Professional output*:
  - Formatted DMS strings (129°33'21.0")
  - Comprehensive tabular computation reports matching field book standards
  - Step-by-step adjustment tracing for auditability
- *Visualization*: Cartesian plot of cut vectors for quality control and geometric validation
- *Zero external dependencies* beyond standard scientific stack (numpy, matplotlib, math)

## 🎯 Use Cases

- Free-station setup in construction surveying
- Traverse network adjustment and verification
- Educational tool for geomatics/surveying coursework
- Field computation validation before instrument setup
- Coordinate determination in GNSS-denied environments

## 📐 Methodology

Implements classical techniques from surveying practice:
- *Intersection geometry* using bearing intersections
- *JOIN calculations* (bearing and distance between points)
- *Polar-to-rectangular transformations* for coordinate propagation
- *Dual cut methodology*:
  - Cut in Northing: ΔN' = cot(bearing) × ΔE
  - Cut in Easting: ΔE' = tan(bearing) × ΔN
- Provisional coordinate averaging from dual intersection solutions

  Designed by,
  kibetcleii807
