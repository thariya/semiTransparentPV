# Coloured and transparent window design using textured semi-transparent solar cell

## Problem Statement
- Current solar cell designs are either opaque or darkly colored
- There is a significant opportunity to install solar panels on building windows and roofs if solar cells could be made to custom colours or tranparent altogether
- However, designing the optics to match that requirements is extremely difficult
- One approach is to interleave transparent glass with strips of solar cell such that the solar cells are invisible to the naked eye
- So as a starting point, can we identify what colours will appear to a human when standing on the other side of a solar panel made with such an interleaved design?

![Image](images/pipeline.png)

## Results
- We use the python meep FDTD library to solve Maxwell's equations for interleaved solar cell stacks with various solar cell widths
- We see that various colours are generated when the sizes of the structures are changed
- Reflected and Transmitted light have different behaviours and different colours

