# GARM-LS source code

This is the source code (CPU version) for our paper: [GARM-LS: A Gradient-Augmented Reference-Map Method for Level-Set Fluid Simulation](https://dl.acm.org/doi/10.1145/3618377)

This code runs on Windows 10 and above. However, it does not contain platform-specific code, so hopefully it can also run on other systems (not tested).

## Usage

The code depends on [xmake](https://xmake.io/#/). Just run `xmake` at the root folder and it will automatically install the required dependencies (Eigen, yaml-cpp, glad, glfw).

To run a test case, the command is as follows:
```
xmake r -w . demo -o output -d <DIMENSION> -t <CASE>
```
For example, running a 2d test case 0 yields the following command:
```
xmake r -w . demo -o output -d 2 -t 0
```
To preview the result, run
```
xmake r -w . viewer -o output
```
Play: `Ctrl+P` Reset: `Ctrl+R`

## Citing

If you use this code please cite:

```
@article{10.1145/3618377,
  author = {Li, Xingqiao and Ni, Xingyu and Zhu, Bo and Wang, Bin and Chen, Baoquan},
  title = {GARM-LS: A Gradient-Augmented Reference-Map Method for Level-Set Fluid Simulation},
  year = {2023},
  volume = {42},
  number = {6},
  journal = {ACM Trans. Graph.},
}
```

## License

MIT License
