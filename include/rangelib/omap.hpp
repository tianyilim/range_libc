#ifndef RANGELIB_OMAP_HPP_
#define RANGELIB_OMAP_HPP_

#include <vector>

#include "RangeUtils.h"
#include "distance_transform.h"
#include "lodepng/lodepng.h"

/// @file This defines the Occupancy Grid class.

namespace ranges {

/// @brief Occupancy Grid
class OMap {
   public:
    using Grid_t = std::vector<std::vector<bool>>;
    using RawGrid_t = std::vector<std::vector<float>>;

    /// @brief Constructor from sizes
    /// @param[in] w width of occupancy grid
    /// @param[in] h height of occupancy grid
    OMap(unsigned w = 0, unsigned h = 0) : _hasError{false}, _width{w}, _height{h}, _filename{""}
    {
        // Initializes a blank occupancy grid map.
        for (int i = 0; i < w; ++i) {
            std::vector<bool> y_axis;
            for (int q = 0; q < h; ++q)
                y_axis.push_back(false);
            _OccupancyGrid.push_back(y_axis);
        }
    }

    /// @brief Constructor from an image file and belief threshold
    /// @param[in] filename path to the image file
    /// @param[in] threshold
    OMap(std::string filename, float threshold = 128) : _hasError(false), _filename(filename)
    {
        unsigned error;
        unsigned char *image;

        error = lodepng_decode32_file(&image, &_width, &_height, filename.c_str());
        if (error) {
            std::cerr << "ERROR " << error << lodepng_error_text(error) << std::endl;
            _hasError = true;
            return;
        }

        for (int i = 0; i < _width; ++i) {
            std::vector<bool> y_axis;
            for (int q = 0; q < _height; ++q)
                y_axis.push_back(false);
            _OccupancyGrid.push_back(y_axis);
        }

        for (int i = 0; i < _width; ++i) {
            std::vector<float> y_axis;
            for (int q = 0; q < _height; ++q)
                y_axis.push_back(0);
            _SignedDistFunc.push_back(y_axis);
        }

        for (int y = 0; y < _height; ++y) {
            for (int x = 0; x < _width; ++x) {
                unsigned idx = 4 * y * _width + 4 * x;
                int r = image[idx + 2];
                int g = image[idx + 1];
                int b = image[idx + 0];
                int gray = (int)utils::rgb2gray(r, g, b);
                if (gray < threshold)
                    _OccupancyGrid[x][y] = true;
                _SignedDistFunc[x][y] = gray;
            }
        }
    }

    /// @brief Return the occupancy at x, y
    inline bool isOccupied(int x, int y)
    {
        if (x < 0 || x >= _width || y < 0 || y >= _height)
            return false;
        else
            return _OccupancyGrid[x][y];
    }

    // TODO this should be moved
    /// @brief Returns the signed distance value at x, y
    inline float signedDistanceValue(int x, int y)
    {
        if (x < 0 || x >= _width || y < 0 || y >= _height)
            return -1;  // this should throw
        else
            return _SignedDistFunc[x][y];
    }

    /// @brief Save map to an image file
    /// @param[in] filename
    /// @return error code
    bool save(std::string &filename)
    {
        std::vector<unsigned char> png;
        lodepng::State state;  // optionally customize this one
        char image[_width * _height * 4];

        for (int y = 0; y < _height; ++y) {
            for (int x = 0; x < _width; ++x) {
                unsigned idx = 4 * y * _width + 4 * x;

                image[idx + 2] = (char)255;
                image[idx + 1] = (char)255;
                image[idx + 0] = (char)255;
                image[idx + 3] = (char)255;

                if (_OccupancyGrid[x][y]) {
                    image[idx + 0] = 0;
                    image[idx + 1] = 0;
                    image[idx + 2] = 0;
                }
            }
        }
        unsigned error = lodepng::encode(png, reinterpret_cast<const unsigned char *>(image), _width, _height, state);
        if (!error)
            lodepng::save_file(png, filename);
        else
            std::cerr << "Image encoder error " << error << ": " << lodepng_error_text(error) << std::endl;

        return error;
    }

    /// @brief Creates an edge map from the existing OMap
    /// @param[in] count_corners
    /// @return
    OMap makeEdgeMap(bool count_corners)
    {
        OMap edge_map = OMap(_width, _height);
        for (int x = 0; x < _width; ++x) {
            for (int y = 0; y < _height; ++y) {
                if (!isOccupied(x, y))
                    continue;

                std::vector<std::pair<int, int>> outline = utils::outline(x, y, count_corners);
                for (int i = 0; i < outline.size(); ++i) {
                    int cx, cy;

                    std::tie(cx, cy) = outline[i];
                    if (0 <= cx && 0 <= cy && cx < _width && cy < _height && !isOccupied(cx, cy)) {
                        edge_map.grid()[x][y] = true;
                        break;
                    }
                }
            }
        }
        return edge_map;
    }

    /// @brief returns memory usage in bytes
    inline int memory()
    {
        return sizeof(bool) * _width * _height;
    }

    /// @brief error accessor
    inline bool error()
    {
        return _hasError;
    }

    ///@brief modifier view into _OccupancyGrid
    Grid_t &grid() { return _OccupancyGrid; }

    /// @brief Const view into _SignedDistFunc
    const RawGrid_t &rawGrid() const { return _SignedDistFunc; }
    ///@brief modifier view into _SignedDistFunc
    RawGrid_t &rawGrid() { return _SignedDistFunc; }

    const unsigned width() const { return _width; }
    void setWidth(unsigned w) { _width = w; }

    const unsigned height() const { return _height; }
    void setHeight(unsigned h) { _height = h; }

    /// @brief filename accessor
    const std::string &filename() const { return _filename; }

    // TODO these should be added to some sort of struct with getter/setter
    float _worldScale;     ///< Real-world values
    float _worldAngle;     ///< Real-world values
    float _worldOriginX;   ///< Real-world values
    float _worldOriginY;   ///< Real-world values
    float _worldSinAngle;  ///< Real-world values
    float _worldCosAngle;  ///< Real-world values

   protected:
    // Members
    bool _hasError;  ///< error flag
    // TODO this should be removed
    unsigned _width;        ///< x axis
    unsigned _height;       ///< y axis
    std::string _filename;  ///< filename of image

    Grid_t _OccupancyGrid;      ///< Grid of boolean occupied/not-occupied values
    RawGrid_t _SignedDistFunc;  ///< Signed Distance Function
    // TODO This should be moved to the DistanceTransform subclass as it is not used
};

class DistanceTransform : public OMap {
   public:
    /// @brief Forwards to OMap constructor
    /// @param[in] w Occupancy grid width
    /// @param[in] h Occupancy grid height
    DistanceTransform(unsigned w = 0, unsigned h = 0) : OMap(w, h) {}

    ///@brief computes the distance transform of a given OMap
    ///@param[in] map Input OMap
    DistanceTransform(OMap &map) : OMap(map)
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
        _SignedDistFunc.clear();
        for (int x = 0; x < _width; ++x) {
            std::vector<float> y_axis;
            for (int y = 0; y < _height; ++y) {
                y_axis.push_back(f[x][y]);
            }
            _SignedDistFunc.push_back(y_axis);
        }
    }
};

}  // namespace ranges

#endif