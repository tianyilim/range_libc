#ifndef RANGELIB_RAY_MARCHING_HPP_
#define RANGELIB_RAY_MARCHING_HPP_

#include "rangelib/omap.hpp"
#include "rangelib/range_method.hpp"

namespace ranges {

class RayMarching : public RangeMethod {
   public:
    RayMarching(DistanceTransform map, float maxRange, float distThreshold = 0.0,
                float stepCoeff = 0.999)
        : RangeMethod(map, maxRange), _distThreshold{distThreshold}, _stepCoeff{stepCoeff}
    {
    }

    float calc_range(const float x, const float y, const float heading) const override
    {
        // Catch illegal array access
        if (x < 0 || x >= _distTransform.width() || y < 0 || y >= _distTransform.height()) {
            return _maxRange;
        }

        float ray_direction_x = cosf(heading);
        float ray_direction_y = sinf(heading);

        float vecLength = 0.0;
        while (vecLength < _maxRange) {
            int px = std::round(x + ray_direction_x * vecLength);
            int py = std::round(y + ray_direction_y * vecLength);

            if (px >= (int)_distTransform.width() || px < 0 || py < 0 ||
                py >= (int)_distTransform.height()) {
                return _maxRange;
            }

            float distToNearestObstacle = _distTransform.signedDistanceValue(px, py);
            if (distToNearestObstacle <= _distThreshold) {
                return std::hypot(px - x, py - y);
            }

            vecLength += std::max<float>(distToNearestObstacle * _stepCoeff, 1.0);
        }
        return _maxRange;
    }

    int memory() const override { return _distTransform.memory(); }

   protected:
    const float _distThreshold;
    const float _stepCoeff;
};

}  // namespace ranges

#endif