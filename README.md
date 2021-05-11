# GCMQ --- OpenSees Implementation

This repository contains the OpenSees implementation of GCMQ element.

The provided `OpenSees.exe` under `binary` folder is compiled based on official source code of version `3.2.2`. Dynamic linked library is not provided due to all kinds of compatibility issues.

Please use the following command to create a new GCMQ element:

```tcl
element GCMQ (1) (2...5) (6) (7) [8]
# (1) int, unique element tag
# (2...5) int, four node tags that define element geometry
# (6) double, thickness of element
# (7) int, tag of material model used
# [8] int, integration scheme tag, 1<=>I, 2<=>L, 3<=>G, default: 1
```

For better functionality, GCMQ is also implemented in [suanPan](https://github.com/TLCFEM/suanPan). Readers are recommended to try it out!
