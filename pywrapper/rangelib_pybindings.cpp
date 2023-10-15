#include <pybind11/cast.h>
#include <pybind11/detail/common.h>
#include <pybind11/eigen.h>
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

// Bindings for base virtual class
class PyRangeMethodBaseClass : public RangeMethod {
   public:
    /* Inherit the constructors */
    using RangeMethod::RangeMethod;

    /* Trampoline code for redirection (need one for each virtual function) */
    float calc_range(float x, float y, float heading) const override
    {
        PYBIND11_OVERRIDE_PURE(float,        /* Return type */
                               RangeMethod,  /* Parent class */
                               calc_range,   /* Name of function in C++ (must match Python name) */
                               x, y, heading /* Argument(s) */
        );
    }

    // Trailing comma needed for macro resolution
    int memory() const override { PYBIND11_OVERRIDE_PURE(int, RangeMethod, memory, ); }
};

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
    py::class_<RangeMethod, PyRangeMethodBaseClass> PyRangeMethod(m, "_RangeMethod");
    PyRangeMethod.def(py::init<DistanceTransform, float>(), py::arg("map"), py::arg("maxRange"),
                      "Range Method pure virtual base class. Do not call this directly.");
    PyRangeMethod.def_property("distTransform", &RangeMethod::getMap, nullptr);
    PyRangeMethod.def_property("maxRange", &RangeMethod::maxRange, nullptr);
    PyRangeMethod.def("batchCalcRange", &RangeMethod::batchCalcRange, py::arg("ins"),
                      py::return_value_policy::move,
                      "Wrapper function to call calc_range repeatedly. Input a 3xn "
                      "numpy array for the inputs and get a 1xn numpy array of the outputs.");
    PyRangeMethod.def("numpy_calc_range_angles", &RangeMethod::numpy_calc_range_angles,
                      py::arg("ins"), py::arg("angles"), py::arg("outs"), py::arg("num_particles"),
                      py::arg("num_angles"),
                      "Wrapper function to call calc_range repeatedly. Indexing assumes 3 x "
                      "num_particles array of sensor poses, 1 x num_angles array of query angles, "
                      "1 x num_angles*num_particles array of output values.");
    PyRangeMethod.def(
        "set_sensor_model", &RangeMethod::set_sensor_model, py::arg("table"),
        py::arg("table_width"),
        "Sets a pre-calculated sensor model (table_width x table_width array). The "
        "sensor model is the likelihood score of observing a particular measurement.");
    PyRangeMethod.def(
        "eval_sensor_model", &RangeMethod::eval_sensor_model, py::arg("observations"),
        py::arg("ranges"), py::arg("outs"), py::arg("rays_per_particle"), py::arg("num_particles"),
        "Evaluating the (discretized) sensor model is equivalent to a 2D array lookup.");
    PyRangeMethod.def("calc_range_repeat_angles_eval_sensor_model",
                      &RangeMethod::calc_range_repeat_angles_eval_sensor_model, py::arg("ins"),
                      py::arg("angles"), py::arg("observations"), py::arg("weights"),
                      py::arg("num_particles"), py::arg("num_angles"),
                      "Function that calculates expected ranges for each particle and directly "
                      "outputs particle weights");
    PyRangeMethod.def("mapCoordinates", &RangeMethod::mapCoordinates, py::arg("x_world"),
                      py::arg("y_world"), py::arg("z_world"),
                      "Convert world coordinates into map-discretized values");

    /// Bresenham's Implementation
    py::class_<BresenhamsLine> PyBresenham(m, "Bresenham", PyRangeMethod);
    PyBresenham.def(py::init<DistanceTransform, float>(), py::arg("map"), py::arg("maxRange"),
                    "Bresenham's Line method constructor");

    /// Ray Marching Implementation
    py::class_<RayMarching> PyRayMarching(m, "RayMarching", PyRangeMethod);
    PyRayMarching.def(py::init<DistanceTransform, float>(), py::arg("map"), py::arg("maxRange"),
                      "Ray Casting method constructor");

    /// Giant Lookup Table implementation
    py::class_<GiantLUTCast> PyLookupTable(m, "GiantLUT", PyRangeMethod);
    PyLookupTable.def(py::init<DistanceTransform, float, unsigned>(), py::arg("map"),
                      py::arg("maxRange"), py::arg("theta_discretization"),
                      "Giant Lookup Table method constructor");
}

}  // namespace ranges