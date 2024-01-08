## Convex Geometric Motion Planning on Lie Groups via Sparse Moment Relaxation

This project solve the motion planning problem using full rigid body dynamics on matrix Lie group. 

To run the code, please install the TSSOS in Julia for global optimization and YALMIP in MATLAB for local search. 
MOSEK is required to solve the semidefinite programming. 

## Requirement:

- [`TSSOS'](https://github.com/wangjie212/TSSOS.git)
- [`MOSEK'](https://www.mosek.com/downloads/)
- [`YALMIP'](https://yalmip.github.io/download/)


## Citation:
This project is released to reproduce the result of the following conference paper and its journal extension. 

```
@INPROCEEDINGS{Teng-RSS-23, 
    AUTHOR    = {Sangli Teng AND Ashkan Jasour AND Ram Vasudevan AND Maani Ghaffari Jadidi}, 
    TITLE     = {{Convex Geometric Motion Planning on Lie Groups via Moment Relaxation}}, 
    BOOKTITLE = {Proceedings of Robotics: Science and Systems}, 
    YEAR      = {2023}, 
    ADDRESS   = {Daegu, Republic of Korea}, 
    MONTH     = {July}, 
    DOI       = {10.15607/RSS.2023.XIX.058} 
} 
```
