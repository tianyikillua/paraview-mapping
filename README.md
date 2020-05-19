# ParaView plugin for mapping finite element data between meshes

This repository contains a ParaView Python plugin (`MappingFilter.py`) that can be loaded by ParaView to transfer finite element data defined on nodes (P1 fields) or on cells (P0 fields) between two non-matching meshes. The actual mapping algorithm is performed by the [MEDCoupling](https://docs.salome-platform.org/latest/dev/MEDCoupling/developer/index.html) library. As a ParaView filter, this plugin provides a handy graphical interface that allows you to
- Pre-process the source and target meshes as well as the field to be transferred from the source mesh
- Perform the mapping of a scalar / vector / tensor field within the GUI
- Post-process the mapped field on the target mesh and export it to an external file

![](https://user-images.githubusercontent.com/4027283/69738195-f6d03300-1135-11ea-8832-5af1950f9bb0.png)

## Installation and updates

If you have downloaded a binary version of ParaView, you may proceed as follows
1. Download the [`medcoupling`](https://github.com/tianyikillua/medcoupling) Python package library and put the `medcoupling` folder into the `site-packages` directory of ParaView. For instance, under Windows, it is `bin\Lib\site-packages`. You need to make sure that ParaView uses a Python version that supports `medcoupling`, that is at least Python 3.7.
2. Download `MappingFilter.py` and load the plugin under ParaView, via *Tools* / *Manage Plugins* / *Load New*. You can optionally check the option *Auto Load*.

To ensure that the current plugin is up to date, you may clone this repository and frequently pull the latest updates
``` sh
git clone https://github.com/tianyikillua/paraview-mapping.git
git pull
```

## Usage

The following animation illustrates the mapping of a scalar field between two 2d meshes. The source mesh is mainly comprised of quadrilateral elements while the target mesh is coarser and is composed of triangular elements.

The main steps are
- Load the source and target meshes in ParaView
- Select the source mesh and load the **Mapping** filter defined by this plugin
- Define and re-check the source and target meshes
- Define mapping parameters. You may refer to the [pymapping](https://github.com/tianyikillua/pymapping) and [MEDCoupling](https://docs.salome-platform.org/latest/dev/MEDCoupling/developer/index.html) libraries for more information.
- If everything goes well, the mapped field is now available in the target mesh

![](https://user-images.githubusercontent.com/4027283/69740883-66482180-113a-11ea-8670-8ec0ec5625da.gif)


## License

`MappingFilter.py` is published under the [MIT license](https://en.wikipedia.org/wiki/MIT_License).
