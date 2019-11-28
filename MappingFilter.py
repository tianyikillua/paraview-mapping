import medcoupling as mc
import numpy as np

from paraview.util.vtkAlgorithm import (
    VTKPythonAlgorithmBase,
    smdomain,
    smproperty,
    smproxy,
)
from vtkmodules.numpy_interface import dataset_adapter as dsa
from vtkmodules.vtkCommonDataModel import vtkUnstructuredGrid

__author__ = "Tianyi Li"
__email__ = "tianyikillua@gmail.com"
__copyright__ = "Copyright (c) 2019 {} <{}>".format(__author__, __email__)
__license__ = "License :: OSI Approved :: MIT License"
__version__ = "0.1.0"
__status__ = "Development Status :: 4 - Beta"

available_fields = [""]

vtk_to_mc_type = {
    1: mc.NORM_POINT1,  # vertex
    3: mc.NORM_SEG2,  # line
    5: mc.NORM_TRI3,  # triangle
    9: mc.NORM_QUAD4,  # quad
    10: mc.NORM_TETRA4,  # tetra
    12: mc.NORM_HEXA8,  # hexahedron
    14: mc.NORM_PYRA5,  # pyramid
}
mc_to_vtk_type = {v: k for k, v in vtk_to_mc_type.items()}

celltype_3d = [10, 12, 14]
celltype_2d = [5, 9]
celltype_1d = [3]
celltype_0d = [1]


def meshdim(cell_types):
    """
    Determine the mesh dimension of a VTK mesh
    """
    cell_types_unique = np.unique(cell_types)
    if len(set(celltype_3d).intersection(cell_types_unique)) > 0:
        meshdim = 3
    elif len(set(celltype_2d).intersection(cell_types_unique)) > 0:
        meshdim = 2
    else:
        meshdim = 1
    return meshdim


def read_point_cell_data(mesh):
    """
    Read point and cell data from a VTK mesh
    """
    # Adapted from https://github.com/nschloe/meshio/blob/master/test/legacy_reader.py
    def _read_data(data):
        out = {}
        for i in range(data.VTKObject.GetNumberOfArrays()):
            name = data.VTKObject.GetArrayName(i)
            array = np.asarray(data.GetArray(i))
            out[name] = array
        return out

    point_data = _read_data(mesh.GetPointData())
    cell_data = _read_data(mesh.GetCellData())
    return point_data, cell_data


def mesh_mc_from_VTK(mesh, check=False):
    """
    Convert a VTK mesh to a MEDCoupling mesh

    TODO: Prune lower-dimensional cells
    """
    # Read points and cells from VTK mesh
    points = np.asarray(mesh.GetPoints()).copy()
    cell_conn = np.asarray(mesh.GetCells()).copy()
    cell_offsets = np.asarray(mesh.GetCellLocations()).copy()
    cell_types = np.asarray(mesh.GetCellTypes()).copy()

    # Initialization
    mesh_mc = mc.MEDCouplingUMesh("mesh", meshdim(cell_types))

    # Points
    if np.isclose(np.linalg.norm(points[:, 2]), 0):
        points = points[:, [0, 1]].copy()
    coords = mc.DataArrayDouble(points)
    mesh_mc.setCoords(coords)

    # Cells
    cell_conn[cell_offsets] = [vtk_to_mc_type[cell_type] for cell_type in cell_types]
    conn = mc.DataArrayInt(cell_conn.astype(np.int32))
    cell_offsets = np.hstack([cell_offsets, [len(cell_conn)]])
    conn_index = mc.DataArrayInt(cell_offsets.astype(np.int32))
    mesh_mc.setConnectivity(conn, conn_index)

    if check:
        mesh_mc.checkConsistency()
    return mesh_mc


def field_mc_from_VTK(
    mesh, field_name, on="points", mesh_mc=None, nature="IntensiveMaximum"
):
    """
    Convert a VTK field to a MEDCoupling field
    """
    assert on in ["points", "cells"]
    if on == "points":
        field = mc.MEDCouplingFieldDouble(mc.ON_NODES, mc.NO_TIME)
    else:
        field = mc.MEDCouplingFieldDouble(mc.ON_CELLS, mc.NO_TIME)
    field.setName(field_name)
    if mesh_mc is None:
        mesh_mc = mesh_mc_from_VTK(mesh)
    field.setMesh(mesh_mc)

    # Read point and cell data
    point_data, cell_data = read_point_cell_data(mesh)
    if on == "points":
        assert field_name in point_data
        field.setArray(mc.DataArrayDouble(point_data[field_name]))
    else:
        assert field_name in cell_data
        field.setArray(mc.DataArrayDouble(cell_data[field_name]))

    field.setNature(eval("mc." + nature))
    return field


@smproxy.filter(label="Mapping", support_reload=False)
@smproperty.input(name="TargetMesh", port_index=1)
@smdomain.datatype(dataTypes=["vtkUnstructuredGrid"], composite_data_supported=False)
@smproperty.input(name="SourceMesh", port_index=0)
@smdomain.datatype(dataTypes=["vtkUnstructuredGrid"], composite_data_supported=False)
class MappingFilter(VTKPythonAlgorithmBase):
    def __init__(self):
        VTKPythonAlgorithmBase.__init__(
            self, nInputPorts=2, nOutputPorts=1, outputType="vtkUnstructuredGrid"
        )
        self._method = "P1P0"
        self._field_name = None
        self._default_value = np.nan
        self._intersection_type = "Automatic"
        self._nature = "IntensiveMaximum"

    # Field to transfer
    @smproperty.stringvector(name="Fields", information_only="1")
    def GetFields(self):
        return available_fields

    @smproperty.stringvector(name="FieldToTransfer", number_of_elements="1")
    @smdomain.xml(
        """
        <StringListDomain name="list">
            <RequiredProperties>
                <Property name="Fields" function="StringInfo"/>
            </RequiredProperties>
        </StringListDomain>
        """
    )
    def SetFieldName(self, field_name):
        if field_name == "":
            field_name = None

        if self._field_name != field_name:
            self._field_name = field_name
            self.Modified()

    # Mapping method
    @smproperty.stringvector(name="Methods", information_only="1")
    def GetMethods(self):
        return ["P1P0", "P1P1", "P0P0", "P0P1"]

    @smproperty.stringvector(name="MappingMethods", number_of_elements="1")
    @smdomain.xml(
        """
        <StringListDomain name="list">
            <RequiredProperties>
                <Property name="Methods" function="StringInfo"/>
            </RequiredProperties>
        </StringListDomain>
        """
    )
    def SetMethod(self, method):
        if self._method != method:
            self._method = method
            self.Modified()

    # Default mapped value
    @smproperty.doublevector(name="DefaultValue", default_values=-1e24)
    @smdomain.doublerange(min=-1e24, max=1e24)
    def SetDefaultValue(self, default_value):
        if np.isclose(default_value, -1e24):
            default_value = np.nan
        if np.isnan(self._default_value) and np.isnan(default_value):
            pass
        elif not np.isclose(self._default_value, default_value):
            self._default_value = default_value
            self.Modified()

    # Intersection algorithm
    @smproperty.stringvector(name="IntersectionTypes", information_only="1")
    def GetIntersectionTypes(self):
        return ["Automatic", "PointLocator", "Triangulation"]

    @smproperty.stringvector(name="IntersectionAlgorithm", number_of_elements="1")
    @smdomain.xml(
        """
        <StringListDomain name="list">
            <RequiredProperties>
                <Property name="IntersectionTypes" function="StringInfo"/>
            </RequiredProperties>
        </StringListDomain>
        """
    )
    def SetIntersectionType(self, intersection_type):
        if self._intersection_type != intersection_type:
            self._intersection_type = intersection_type
            self.Modified()

    # Nature of the field to transfer
    @smproperty.stringvector(name="Natures", information_only="1")
    def GetNatures(self):
        return [
            "IntensiveMaximum",
            "IntensiveConservation",
            "ExtensiveMaximum",
            "ExtensiveConservation",
        ]

    @smproperty.stringvector(name="NatureOfTheField", number_of_elements="1")
    @smdomain.xml(
        """
        <StringListDomain name="list">
            <RequiredProperties>
                <Property name="Natures" function="StringInfo"/>
            </RequiredProperties>
        </StringListDomain>
        """
    )
    def SetNature(self, nature):
        if self._nature != nature:
            self._nature = nature
            self.Modified()

    def RequestData(self, request, inInfo, outInfo):
        global available_fields

        mesh_source = dsa.WrapDataObject(vtkUnstructuredGrid.GetData(inInfo[0]))
        mesh_target = dsa.WrapDataObject(vtkUnstructuredGrid.GetData(inInfo[1]))
        output = dsa.WrapDataObject(vtkUnstructuredGrid.GetData(outInfo))
        output.ShallowCopy(vtkUnstructuredGrid.GetData(inInfo[1]))

        # Update available fields to transfer from the source mesh
        point_data, cell_data = read_point_cell_data(mesh_source)
        available_fields = list(set(point_data.keys()).union(set(cell_data.keys())))
        self.Modified()

        # If only one available field exists, choose this field by default
        if len(available_fields) == 1:
            self._field_name = available_fields[0]

        # Perform mapping
        if self._field_name is not None:
            mapper = mc.MEDCouplingRemapper()

            # Construct MEDCoupling meshes
            mesh_source_mc = mesh_mc_from_VTK(mesh_source)
            mesh_target_mc = mesh_mc_from_VTK(mesh_target)

            # Preparation
            if self._intersection_type != "Automatic":
                mapper.setIntersectionType(eval("mc." + self._intersection_type))
            elif self._method[:2] == "P1":
                mapper.setIntersectionType(mc.PointLocator)
            mapper.prepare(mesh_source_mc, mesh_target_mc, self._method)

            # Construct field to transfer
            if self._method[:2] == "P1":
                on = "points"
            else:
                on = "cells"
            field_source = field_mc_from_VTK(
                mesh_source,
                self._field_name,
                on=on,
                mesh_mc=mesh_source_mc,
                nature=self._nature,
            )

            # Mapping
            field_target = mapper.transferField(
                field_source, dftValue=self._default_value
            )
            array = field_target.getArray().toNumPyArray()

            # Export mapped field to output
            if self._method[2:] == "P1":
                output.PointData.append(array, self._field_name)
            else:
                output.CellData.append(array, self._field_name)

        return 1
