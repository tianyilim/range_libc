#ifndef RANGELIB_GIANT_LOOKUP_TABLE_HPP_
#define RANGELIB_GIANT_LOOKUP_TABLE_HPP_

#include "rangelib/omap.hpp"
#include "rangelib/range_method.hpp"

#define _GIANT_LUT_SHORT_DATATYPE 1

namespace ranges {
class GiantLUTCast : public RangeMethod {
   public:
#if _GIANT_LUT_SHORT_DATATYPE
    typedef uint16_t lut_t;
#else
    typedef float lut_t;
#endif

    /// @brief Returns a LUT index based on the discretized inputs
    inline size_t getLutIdx(size_t x, size_t y, size_t theta)
    {
        return theta + y * _thetaDiscretization + x * _thetaDiscretization * _height;
    }

    /// @brief Constructor.
    /// @param[in] m Input Occupancy Grid Map
    /// @param[in] mr Max range
    /// @param[in] td theta discretization
    GiantLUTCast(OMap m, float mr, int td)
        : _height{m.height()},
          _width{m.width()},
          _thetaDiscretization{td},
          RangeMethod(m, mr),
          _thetaDiscretization_div_M_2PI{(float)_thetaDiscretization / M_2PI},
          _M_2PI_div_thetaDiscretization{M_2PI / (float)_thetaDiscretization},
          _maxDivLimits{max_range / std::numeric_limits<uint16_t>::max()},
          _limitsDivMax{std::numeric_limits<uint16_t>::max() / max_range}
    {
        // To initalize the LUT, use Ray Marching
        RayMarching seed_cast = RayMarching(m, mr);

        _GiantLUT.resize(_height * _width * _thetaDiscretization);
        for (int x = 0; x < _width; ++x) {
            for (int y = 0; y < _height; ++y) {
                for (int i = 0; i < _thetaDiscretization; ++i) {
                    float angle = i * _M_2PI_div_thetaDiscretization;
                    float r = seed_cast.calc_range(x, y, angle);

#if _GIANT_LUT_SHORT_DATATYPE
                    r = std::min(max_range, r);
                    uint16_t val = r * _limitsDivMax;
                    _GiantLUT[getLutIdx(x, y, i)] = val;
#else
                    _GiantLUT[getLutIdx(x, y, i)] = r;
#endif
                }
            }
        }
    }

    int lut_size() { return map.width() * map.height() * _thetaDiscretization * sizeof(lut_t); }

    int memory() { return lut_size(); }

    /// @brief takes a continuous theta space and returns the nearest theta in the discrete LUT
    /// space as well as the bin index that the given theta falls into
    /// @param[in] theta the angle
    /// @return the discretized value
    int discretize_theta(float theta)
    {
        theta = fmod(theta, M_2PI);
        // fmod does not wrap the angle into the positive range, so this will fix that if necessary
        if (theta < 0.0) theta += M_2PI;

#if _USE_FAST_ROUND == 1
        int rounded = int(theta * _thetaDiscretization_div_M_2PI + 0.5);
#else
        int rounded = (int)roundf(theta * _thetaDiscretization_div_M_2PI);
#endif
        int binned = rounded % _thetaDiscretization;

        return binned;
    }

    float calc_range(float x, float y, float heading)
    {
        if (x < 0 || x >= map.width() || y < 0 || y >= map.height()) {
            return max_range;
        }

        size_t idx = getLutIdx((size_t)x, (size_t)y, discretize_theta(heading));
#if _GIANT_LUT_SHORT_DATATYPE
        return _GiantLUT[idx] * _maxDivLimits;
#else
        return _GiantLUT[idx];
#endif
    }

    /// @brief Returns a slice of the LUT based on a theta value
    /// @param[in] theta the slice value
    /// @return a DistanceTransform unique_ptr
    std::unique_ptr<DistanceTransform> get_slice(float theta)
    {
        std::unique_ptr<DistanceTransform> slice =
            std::make_unique<DistanceTransform>(_width, _height);
        const int dtheta = discretize_theta(theta);

        size_t idx = dtheta;  // iterator index

        for (int x = 0; x < _width; ++x) {
            for (int y = 0; y < _height; ++y) {
                slice->grid()[x][y] = _GiantLUT[idx];
                idx += _thetaDiscretization;
            }
            idx += _thetaDiscretization * _height;
        }

        return slice;
    }

   protected:
    const int _height, _width;                   ///< dimensions of omap
    const int _thetaDiscretization;              ///< how many values of theta are inside
    const float _thetaDiscretization_div_M_2PI;  ///< theta_discretization / 2pi
    const float _M_2PI_div_thetaDiscretization;  ///< 2pi / theta_discretization
    const float _maxDivLimits;                   ///< Max range / uint16_t max
    const float _limitsDivMax;                   ///< uint16_t max / Max range

    std::vector<lut_t> _GiantLUT;  ///< actual LUT
};
}  // namespace ranges

#endif