# GCMQ --- OpenSees Implementation

This repository contains the OpenSees implementation of GCMQ element.

The element is available as dynamic link library (x64 only) that is callable by the OpenSees main program. The library file and dependent OpenBLAS library can be found in `/binary` folder. The library is compiled in accordance with OpenSees 2.5.0, if it does not work with other versions please let me know by creating issues.

To use the library, please download and place `GCMQ.dll` and `libopenblas.dll` in the folder contains `OpenSees.exe`. Please use the following command to create a new GCMQ element:

```tcl
element GCMQ (1) (2...5) (6) (7) [8]
# (1) int, unique element tag
# (2...5) int, four node tags that define element geometry
# (6) double, thickness of element
# (7) int, tag of material model used
# [8] int, integration scheme tag, 1<=>I, 2<=>L, 3<=>G, default: 1
```

