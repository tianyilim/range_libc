#ifndef RANGELIB_OMAP_HPP_
#define RANGELIB_OMAP_HPP_

#include <vector>

#include "RangeUtils.h"
#include "distance_transform.h"
#include "lodepng/lodepng.h"

/// @file This defines the Occupancy Grid class.

namespace ranges {

struct WorldValues {
    /// @brief Value constructor
    WorldValues(float worldScale_ = 1.0, float worldAngle_ = 0.0, float worldOriginX_ = 0.0,
                float worldOriginY_ = 0.0)
        : worldScale{worldScale_},
          worldAngle{worldAngle_},
          worldOriginX{worldOriginX_},
          worldOriginY{worldOriginY_},
          worldSinAngle{sinf(worldAngle_)},
          worldCosAngle{cosf(worldAngle_)}
    {
    }

    // Implement default values
    const float worldScale;     ///< Resolution of the map, meters / pixel
    const float worldAngle;     ///< Yaw (ccw rotation) of the lower-left pixel in the map
    const float worldOriginX;   ///< X pose of the lower-left pixel in the map
    const float worldOriginY;   ///< Y pose of the lower-left pixel in the map
    const float worldSinAngle;  ///< Sine of worldAngle
    const float worldCosAngle;  ///< Cosine of worldAngle

    inline bool operator==(const WorldValues &other) const
    {
        return worldScale == other.worldScale && worldAngle == other.worldAngle &&
               worldOriginX == other.worldOriginX && worldOriginY == other.worldOriginY &&
               worldSinAngle == other.worldSinAngle && worldCosAngle == other.worldCosAngle;
    }

    inline bool operator!=(const WorldValues &other) const { return !(*this == other); }
};

/// @brief Occupancy Grid
class OMap {
   public:
    using Grid_t = std::vector<std::vector<bool>>;
    using DistanceFunction_t = std::vector<std::vector<float>>;

    /// @brief Constructs empty map from sizes
    /// @param[in] w width of occupancy grid
    /// @param[in] h height of occupancy grid
    OMap(unsigned w = 0, unsigned h = 0) : _width{w}, _height{h}, _filename{""}
    {
        // Initializes a blank occupancy grid map.
        for (unsigned i = 0; i < w; ++i) {
            std::vector<bool> y_axis;
            for (unsigned q = 0; q < h; ++q) {
                y_axis.push_back(false);
            }
            _OccupancyGrid.push_back(y_axis);
        }
    }

    /// @brief Constructor from an image file
    /// @param[in] filename path to the image file
    /// @param[in] threshold occupancy threshold
    OMap(const std::string &filename, float threshold = 128) : _filename(filename)
    {
        unsigned error;
        unsigned char *image;

        error = lodepng_decode32_file(&image, &_width, &_height, filename.c_str());
        if (error) {
            std::cerr << "ERROR " << error << lodepng_error_text(error) << std::endl;
            throw std::runtime_error("LodePng Error");
        }

        // Initialize blank occupancy grid map
        for (unsigned i = 0; i < _width; ++i) {
            std::vector<bool> y_axis;
            for (unsigned q = 0; q < _height; ++q) {
                y_axis.push_back(false);
            }
            _OccupancyGrid.push_back(y_axis);
        }

        for (unsigned y = 0; y < _height; ++y) {
            for (unsigned x = 0; x < _width; ++x) {
                unsigned idx = 4 * y * _width + 4 * x;
                int r = image[idx + 2];
                int g = image[idx + 1];
                int b = image[idx + 0];
                int gray = (int)utils::rgb2gray(r, g, b);

                // Threshold determines occupancy
                _OccupancyGrid[x][y] = gray < threshold;
            }
        }
    }

    /// @brief Copy Constructor
    OMap(const OMap &other)
        : _width{other.width()},
          _height{other.height()},
          _filename{other.filename()},
          _OccupancyGrid(other.grid())
    {
        _worldValuesPtr = std::make_unique<WorldValues>(
            other._worldValuesPtr->worldScale, other._worldValuesPtr->worldAngle,
            other._worldValuesPtr->worldOriginX, other._worldValuesPtr->worldOriginY);
    }

    /// @brief Move constructor
    OMap(OMap &&other) : _width{other.width()}, _height{other.height()}, _filename{other.filename()}
    {
        _OccupancyGrid = std::move(other.grid());
        _worldValuesPtr.swap(other._worldValuesPtr);
    }

    /// @brief Copy assignment Operator
    OMap &operator=(const OMap &other)
    {
        if (this != &other)  // not a self-assignment
        {
            _width = other.width();
            _height = other.height();
            _filename = other.filename();

            _OccupancyGrid = other.grid();

            _worldValuesPtr = std::make_unique<WorldValues>(
                other._worldValuesPtr->worldScale, other._worldValuesPtr->worldAngle,
                other._worldValuesPtr->worldOriginX, other._worldValuesPtr->worldOriginY);
        }

        return *this;
    }

    /// @brief Return the occupancy at x, y
    inline bool isOccupied(unsigned x, unsigned y) const
    {
        if (x >= _width || y >= _height)
            throw std::invalid_argument("invalid x/y argument to occupancy grid");
        else
            return _OccupancyGrid[x][y];
    }

    /// @brief Save map to an image file
    /// @param[in] filename
    /// @return If there was an error.
    bool save(const std::string &filename)
    {
        std::vector<unsigned char> png;
        lodepng::State state;  // optionally customize this one

        char *image;
        image = new char[_width * _height * 4];

        for (unsigned y = 0; y < _height; ++y) {
            for (unsigned x = 0; x < _width; ++x) {
                unsigned idx = 4 * y * _width + 4 * x;

                image[idx + 2] = (char)255;
                image[idx + 1] = (char)255;
                image[idx + 0] = (char)255;
                image[idx + 3] = (char)255;

                if (isOccupied(x, y)) {
                    image[idx + 0] = 0;
                    image[idx + 1] = 0;
                    image[idx + 2] = 0;
                }
            }
        }
        unsigned error = lodepng::encode(png, reinterpret_cast<const unsigned char *>(image),
                                         _width, _height, state);
        delete[] image;

        if (!error)
            lodepng::save_file(png, filename);
        else
            std::cerr << "Image encoder error " << error << ": " << lodepng_error_text(error)
                      << std::endl;

        return error;
    }

    /// @brief Creates an edge map from the existing OMap. An Edge Map sets occupancy to be 1 only
    /// at edges.
    /// @param[in] count_corners if corners of candidate points should be included in the edge map
    /// @return Edge Map
    OMap makeEdgeMap(bool count_corners = true)
    {
        OMap edge_map = OMap(_width, _height);

        for (unsigned x = 0; x < _width; ++x) {
            for (unsigned y = 0; y < _height; ++y) {
                if (!isOccupied(x, y)) continue;

                // candidate outlines
                std::vector<std::pair<int, int>> outline = utils::outline(x, y, count_corners);
                for (size_t i = 0; i < outline.size(); ++i) {
                    int cx, cy;

                    std::tie(cx, cy) = outline[i];

                    // A point is an edge if it is within bounds and not
                    if (0 <= cx && 0 <= cy && cx < int(_width) && cy < int(_height) &&
                        !isOccupied(cx, cy)) {
                        edge_map.grid()[x][y] = true;
                        break;
                    }
                }
            }
        }
        return edge_map;
    }

    /// @brief returns memory usage in bytes
    virtual inline int memory() const { return sizeof(bool) * _width * _height; }

    ///@brief accessor view into _OccupancyGrid
    const Grid_t &grid() const { return _OccupancyGrid; }
    ///@brief modifier view into _OccupancyGrid
    Grid_t &grid() { return _OccupancyGrid; }

    const Grid_t &getGrid() { return grid(); }
    void setGrid(const Grid_t &gridIn) { _OccupancyGrid = gridIn; }

    unsigned width() const { return _width; }

    unsigned height() const { return _height; }

    /// @brief filename accessor
    const std::string &filename() const { return _filename; }

    /// @brief World values const view
    const WorldValues &worldValues() const { return *_worldValuesPtr; }
    /// @brief setter function for _worldValues
    void setWorldValues(const WorldValues &wv)
    {
        _worldValuesPtr = std::make_unique<WorldValues>(wv.worldScale, wv.worldAngle,
                                                        wv.worldOriginX, wv.worldOriginY);
    }

   protected:
    // Members
    unsigned _width;        ///< x axis
    unsigned _height;       ///< y axis
    std::string _filename;  ///< filename of image

    Grid_t _OccupancyGrid;  ///< Grid of boolean occupied/not-occupied values
    std::unique_ptr<WorldValues> _worldValuesPtr =
        std::make_unique<WorldValues>();  ///< Real world parameters of the map
};

class DistanceTransform : public OMap {
   public:
    /// @brief Forwards to OMap constructor
    /// @param[in] w Occupancy grid width
    /// @param[in] h Occupancy grid height
    DistanceTransform(unsigned w = 0, unsigned h = 0) : OMap(w, h)
    {
        // ensure that _SDF has the same dimension as OMap
        for (unsigned i = 0; i < w; ++i) {
            std::vector<float> y_axis;
            for (unsigned q = 0; q < h; ++q) {
                y_axis.push_back(0);
            }
            _SDF.push_back(y_axis);
        }
    }

    ///@brief computes the distance transform of a given OMap
    ///@param[in] map Input OMap
    DistanceTransform(const OMap &map) : OMap(map)
    {
        std::vector<std::size_t> grid_size({_width, _height});
        dt::MMArray<float, 2> f(grid_size.data());
        dt::MMArray<std::size_t, 2> indices(grid_size.data());

        for (std::size_t i = 0; i < _width; ++i) {
            for (std::size_t j = 0; j < _height; ++j) {
                if (map.isOccupied(i, j)) {
                    f[i][j] = 0.0f;
                }
                else {
                    f[i][j] = std::numeric_limits<float>::max();
                }
            }
        }

        dt::DistanceTransform::distanceTransformL2(f, f, indices, false);

        // allocate space in the vectors
        _SDF.clear();
        for (unsigned x = 0; x < _width; ++x) {
            std::vector<float> y_axis;
            for (unsigned y = 0; y < _height; ++y) {
                y_axis.push_back(f[x][y]);
            }
            _SDF.push_back(y_axis);
        }
    }

    /// @brief Returns the signed distance value at x, y
    inline float signedDistanceValue(unsigned x, unsigned y) const
    {
        if (x >= _width || y >= _height) {
            throw std::invalid_argument("invalid x/y argument to signed distance function");
        }
        else {
            return _SDF[x][y];
        }
    }

    /// @brief Returns the memory usage of the DistanceTransform object.
    /// @return Sum of Occupancy Grid memory and SDF memory.
    inline int memory() const override
    {
        int omapSize = sizeof(bool) * _OccupancyGrid.size();
        int sdfSize = sizeof(float) * _SDF.size();

        return omapSize + sdfSize;
    }

    /// @brief Const view into _SignedDistFunc
    const DistanceFunction_t &SDF() const { return _SDF; }
    const DistanceFunction_t &getSDF() const { return SDF(); }

    ///@brief modifier view into _SignedDistFunc
    DistanceFunction_t &SDF() { return _SDF; }

   private:
    DistanceFunction_t _SDF;
};

}  // namespace ranges

#endif