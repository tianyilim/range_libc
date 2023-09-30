#ifndef RANGELIB_BRESENHAMS_HPP_
#define RANGELIB_BRESENHAMS_HPP_

#include "rangelib/omap.hpp"
#include "rangelib/range_method.hpp"

namespace ranges {

class BresenhamsLine : public RangeMethod {
   public:
    BresenhamsLine(OMap m, float mr) : RangeMethod(m, mr){};

    float calc_range(float x, float y, float heading)
    {
        // first check if the cell underneath the query point is occupied, if so return
        if (_distTransform.isOccupied((int)x, (int)y)) {
            return 0.0;
        }

        /*
         this defines the coordinate system such that
            ------> +x
            |
            |
            \/
            +y
          0* heading lies along the x axis, positive heading rotates towards the positive y axis
        */
        float x0 = y;
        float y0 = x;
        float x1 = y + _maxRange * sinf(heading);
        float y1 = x + _maxRange * cosf(heading);

        bool steep = false;
        if (std::abs(y1 - y0) > std::abs(x1 - x0)) steep = true;

        if (steep) {
            float tmp = x0;
            x0 = y0;
            y0 = tmp;
            tmp = x1;
            x1 = y1;
            y1 = tmp;
        }

        float deltax = std::abs(x1 - x0);
        float deltay = std::abs(y1 - y0);

        float error = 0;
        float deltaerr = deltay;
        float _x = x0;
        float _y = y0;

        int xstep = -1;
        if (x0 < x1) xstep = 1;

        int ystep = -1;
        if (y0 < y1) ystep = 1;

        unsigned width = _distTransform.width();
        unsigned height = _distTransform.height();

        while ((int)_x != (int)(x1 + xstep)) {
            _x += xstep;
            error += deltaerr;

            if (error * 2.00 >= deltax) {
                _y += ystep;
                error -= deltax;
            }

            if (!steep) {
                if (0 <= _y && _y < width && 0 <= _x && _x < height &&
                    _distTransform.isOccupied(_y, _x)) {
                    float xd = _x - x0;
                    float yd = _y - y0;
                    return sqrtf(xd * xd + yd * yd);
                }
            }
            else {
                if (0 <= _x && _x < width && 0 <= _y && _y < height &&
                    _distTransform.isOccupied(_x, _y)) {
                    float xd = _x - x0;
                    float yd = _y - y0;
                    return sqrtf(xd * xd + yd * yd);
                }
            }
        }
        return _maxRange;
    }

    int memory() const override { return _distTransform.memory(); }
};

}  // namespace ranges

#endif