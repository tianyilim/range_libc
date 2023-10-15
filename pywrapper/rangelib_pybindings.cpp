#include <pybind11/detail/common.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <string>

#include "rangelib/bresenhams.hpp"
#include "rangelib/lookup_table.hpp"
#include "rangelib/omap.hpp"
#include "rangelib/range_method.hpp"
#include "rangelib/ray_casting.hpp"

namespace py = pybind11;

namespace ranges {

// Registers all rangelib C++ functions in Python.
PYBIND11_MODULE(rangelib, m)
{
    m.doc() = "The RangeLib Python bindings for C++ code";

    /// World Values
    py::class_<WorldValues> PyWorldValues(m, "WorldValues");
    PyWorldValues.def(py::init<float, float, float, float>(), py::arg("worldScale") = 1.0,
                      py::arg("worldAngle") = 0.0, py::arg("worldOriginX") = 0.0,
                      py::arg("worldOriginY") = 0.0,
                      "Maps real-world quantities to an occupancy grid representation.");
    PyWorldValues.def_readonly("worldScale", &WorldValues::worldScale,
                               "Resolution of the map, meters / pixel");
    PyWorldValues.def_readonly("worldAngle", &WorldValues::worldAngle,
                               "Yaw (ccw rotation) of the lower-left pixel in the map");
    PyWorldValues.def_readonly("worldOriginX", &WorldValues::worldOriginX,
                               "X pose of the lower-left pixel in the map");
    PyWorldValues.def_readonly("worldOriginY", &WorldValues::worldOriginY,
                               "Y pose of the lower-left pixel in the map");
    PyWorldValues.def("__repr__", [](const WorldValues &a) {
        return "<rangelib.WorldValues with worldScale:" + std::to_string(a.worldScale) +
               " worldAngle: " + std::to_string(a.worldAngle * M_PI / 180.0) +
               " deg worldOriginX: " + std::to_string(a.worldOriginX) +
               " worldOriginY: " + std::to_string(a.worldOriginY) + " >";
    });

    /// OMap Class
    py::class_<OMap> PyOMap(m, "OMap");
    PyOMap.def(py::init<unsigned, unsigned>(), py::arg("width") = 0, py::arg("height") = 0,
               "Construct an OMap with a width and a height.");
    PyOMap.def(py::init<const std::string &, float>(), py::arg("filepath"),
               py::arg("tolerance") = 128.0,
               "Construct an OMap from a PNG filepath and a occupancy tolerance.");
    // PyOMap.def(py::init<const OMap &>(), "Copy construct an OMap.");
    // PyOMap.def(py::init<OMap &&>(), "Move construct an OMap.");
    PyOMap.def("isOccupied", &OMap::isOccupied, py::arg("x"), py::arg("y"),
               "Return the occupancy at x, y");
    PyOMap.def("save", &OMap::save, py::arg("filename"), "Save map to an image file");
    PyOMap.def("makeEdgeMap", &OMap::makeEdgeMap, py::arg("countCorners") = true,
               "Creates an edge map from the existing OMap. An Edge Map sets occupancy to "
               "be 1 only at edges.");
    PyOMap.def_property("memory", &OMap::memory, nullptr, "Memory usage in bytes");
    PyOMap.def_property("occupancyGrid", &OMap::getGrid, &OMap::setGrid,
                        "The underlying distance function (2D array of bool)");
    PyOMap.def_property("width", &OMap::width, nullptr);
    PyOMap.def_property("height", &OMap::height, nullptr);
    PyOMap.def_property("filename", &OMap::filename, nullptr);
    PyOMap.def_property("worldValues", &OMap::worldValues, &OMap::setWorldValues,
                        "Real world vals of the OMap");
    PyOMap.def("__repr__", [](const OMap &a) {
        return "<rangelib.OMap with width: " + std::to_string(a.width()) +
               " height: " + std::to_string(a.height()) + " filename: " + a.filename() + " >";
    });
    PyOMap.doc() = "Occupancy Grid Map from rangelib";

    /// Distance Transform Class
    py::class_<DistanceTransform> PyDistanceTransform(m, "DistanceTransform", PyOMap);
    PyDistanceTransform.def(py::init<unsigned, unsigned>(), py::arg("width") = 0,
                            py::arg("height") = 0,
                            "Construct a Distance Transform with a width and a height.");
    PyDistanceTransform.def(py::init<const OMap &>(), py::arg("OMap"),
                            "Construct a Distance Transform from an OMap.");
    PyDistanceTransform.def("signedDistanceValue", &DistanceTransform::signedDistanceValue,
                            py::arg("x"), py::arg("y"),
                            "Query the signed distance value at an index.");
    PyDistanceTransform.def_property("memory", &DistanceTransform::memory, nullptr,
                                     "The memory usage of the Distance Transform in bytes");
    PyDistanceTransform.def_property("distanceFunction", &DistanceTransform::getSDF, nullptr,
                                     "The underlying distance function (2D array of float)");
    PyDistanceTransform.doc() = "Distance Transform from rangelib";

    /// Range Method class
    py::class_<RangeMethod> PyRangeMethod(m, "_RangeMethod");
    PyRangeMethod.def(py::init<OMap, float>());

    /// Bresenham's Implementation
    py::class_<BresenhamsLine> PyBresenham(m, "Bresenham", PyRangeMethod);

    /// Ray Marching Implementation
    py::class_<RayMarching> PyRayMarching(m, "RayMarching", PyRangeMethod);

    /// Giant Lookup Table implementation
    py::class_<GiantLUTCast> PyLookupTable(m, "GiantLUT", PyRangeMethod);
}

}  // namespace ranges