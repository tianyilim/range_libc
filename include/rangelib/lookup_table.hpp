#ifndef RANGELIB_GIANT_LOOKUP_TABLE_HPP_
#define RANGELIB_GIANT_LOOKUP_TABLE_HPP_

#include <cmath>

#include "rangelib/omap.hpp"
#include "rangelib/range_method.hpp"
#include "rangelib/ray_casting.hpp"

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
    inline size_t getLutIdx(size_t x, size_t y, size_t theta) const
    {
        return theta + y * _thetaDiscretization + x * _thetaDiscretization * _lutHeight;
    }

    /// @brief Constructor.
    /// @param[in] m Input Occupancy Grid Map
    /// @param[in] mr Max range
    /// @param[in] td theta discretization
    GiantLUTCast(DistanceTransform map, float maxRange, unsigned thetaDiscretization)
        : RangeMethod(map, maxRange),
          _lutHeight{map.height()},
          _lutWidth{map.width()},
          _thetaDiscretization{thetaDiscretization},
          _thetaDiscretization_div_M_2PI{float(_thetaDiscretization / (M_PIf * 2))},
          _M_2PI_div_thetaDiscretization{float((M_PIf * 2) / _thetaDiscretization)},
          _maxDivLimits{_maxRange / std::numeric_limits<uint16_t>::max()},
          _limitsDivMax{std::numeric_limits<uint16_t>::max() / _maxRange}
    {
        // To initalize the LUT, use Ray Marching
        RayMarching seed_cast = RayMarching(map, maxRange);

        _GiantLUT.resize(_lutHeight * _lutWidth * _thetaDiscretization);
        for (unsigned x = 0; x < _lutWidth; ++x) {
            for (unsigned y = 0; y < _lutHeight; ++y) {
                for (unsigned i = 0; i < _thetaDiscretization; ++i) {
                    float angle = i * _M_2PI_div_thetaDiscretization;
                    float r = seed_cast.calc_range(x, y, angle);

#if _GIANT_LUT_SHORT_DATATYPE
                    r = std::min(_maxRange, r);
                    uint16_t val = r * _limitsDivMax;
                    _GiantLUT[getLutIdx(x, y, i)] = val;
#else
                    _GiantLUT[getLutIdx(x, y, i)] = r;
#endif
                }
            }
        }
    }

    /// @brief returns the memory usage of the internal LUT.
    int lut_size() const { return _lutHeight * _lutWidth * _thetaDiscretization * sizeof(lut_t); }

    int memory() const override { return lut_size() + _distTransform.memory(); }

    /// @brief takes a continuous theta space and returns the nearest theta in the discrete LUT
    /// space as well as the bin index that the given theta falls into
    /// @param[in] theta the angle
    /// @return the discretized value
    int discretize_theta(float theta) const
    {
        theta = fmod(theta, M_2PI);
        // fmod does not wrap the angle into the positive range, so this will fix that if necessary
        if (theta < 0.0) theta += M_2PI;

#if _USE_FAST_ROUND == 1
        int rounded = int(theta * _thetaDiscretization_div_M_2PI + 0.5);
#else
        int rounded = (int)roundf(theta * _thetaDiscretization_div_M_2PI);
#endif
        // perhaps keep angle between (0, _thetaDiscretization-1).
        return rounded % _thetaDiscretization;
    }

    float calc_range(const float x, const float y, const float heading) const override
    {
        // check if this is still the case with different scaling
        if (x < 0 || x >= _lutWidth || y < 0 || y >= _lutHeight) {
            return _maxRange;
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
    std::unique_ptr<DistanceTransform> get_slice(float theta) const
    {
        std::unique_ptr<DistanceTransform> slice =
            std::make_unique<DistanceTransform>(_lutWidth, _lutHeight);
        const int dtheta = discretize_theta(theta);

        size_t idx = dtheta;  // iterator index

        for (unsigned x = 0; x < _lutWidth; ++x) {
            for (unsigned y = 0; y < _lutHeight; ++y) {
                slice->grid()[x][y] = _GiantLUT[idx];
                idx += _thetaDiscretization;
            }
            idx += _thetaDiscretization * _lutHeight;
        }

        return slice;
    }

    unsigned lutHeight() const { return _lutHeight; }
    unsigned lutWidth() const { return _lutWidth; }
    unsigned lutThetaDiscretization() const { return _thetaDiscretization; }
    unsigned lutArraySize() const { return _GiantLUT.size(); }

   protected:
    const unsigned _lutHeight, _lutWidth;        ///< dimensions of omap
    const unsigned _thetaDiscretization;         ///< how many values of theta are inside
    const float _thetaDiscretization_div_M_2PI;  ///< theta_discretization / 2pi
    const float _M_2PI_div_thetaDiscretization;  ///< 2pi / theta_discretization
    const float _maxDivLimits;                   ///< Max range / uint16_t max
    const float _limitsDivMax;                   ///< uint16_t max / Max range

    std::vector<lut_t> _GiantLUT;  ///< actual LUT
};
}  // namespace ranges

#endif