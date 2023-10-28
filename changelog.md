## TODO
### [RangeUtils](./include/RangeUtils.h)
- [ ] Add docs to each of the methods.

### [OMap](./include/rangelib/omap.hpp)
- [ ] Test the `save` function in some way
- [ ] Test edge Map functionality in some way

### [Range Method](./include/rangelib/range_method.hpp)
- [ ] Test sensor model functions
- [ ] Test `numpy_` functions for each range method
- [ ] Rename `numpy_` to be "batch" functions
- [ ] Put sensor model into 1d vector for faster indexing
- [ ] Change function signatures from vectors to work with pybind code (currently broken)

### [Bresenhams](./include/rangelib/bresenhams.hpp)
- [ ] Add test for Map Coordinate transform

### [Ray Casting](./include/rangelib/ray_casting.hpp)

### [Giant LUT](./include/rangelib/lookup_table.hpp)

### [Python Bindings](./pywrapper/rangelib_pybindings.cpp)
- [ ] (automated) testing for the wrapped python functions

### CUDA Support
- [ ] Clean up Ray Casting GPU
- [ ] Add instructions to compile with CUDA

### Testing
- [ ] Edit CMake such that tests are only built when TESTING flag is set.
- [ ] Edit tests such that they always point to the correct directories to open lodepng

### Functionality:
- [ ] Better timing benchmarks that simulate real-world operation
- [ ] Visualizations of the different sensor models and results they might give
- [ ] Wrap particle filter base class into C++ code

### Housekeeping
- [ ] Move implementation into `cpp` files
- [ ] Rename all headers to `hpp`

## 28/10/2023
- Moved the build artefacts from `setup.py` into the `pywrapper` folder
- Update the interfaces for each method to run with Eigen + Numpy, add docs
- Update tests with these interfaces

## 16/10/2023
- Integrate `pybind11` with a `setup.py` which installs rangelib on the system
- Add instructions to install these fiels

## 15/10/2023
- Fix theta error with ray marching (needed to call round instead of straight int cast)
- Fix GLT test failing by increasing the theta_discretization
- Add Python wrapper using `Pybind11`.

## 08/10/2023
- Add extra tests to the repo
- Conversion from WORLD to MAP frame was overhauled. Not sure why / how the previous version worked
- Add analytic calculation for specific scenarios
  - Can use this to quantify the noise of pixelation, and in the future, localization in general.

## 01/10/2023
- Added a copy assignment constructor for OMap and made WorldValues members const.
- Will not add copy assignment to `RangeMethod` because otherwise const members cannot be made const.
- Change Bresenham's implementation such that `calc_range` returns max range (like other methods)

## 30/09/2023
- Added tests.
### Testing Instructions
1. From directory root: `cmake -S . -B build` This configures the CMakeLists. But, VS Code also does this automatically on saving CMakeLists.txt.
2. From the `build` folder: `cmake --build .` This builds all targets, which will then be put into `build/bin`
3. From `build/bin`: Run the relevant test targets.

## 24/09/2023
### Changes from original `rangelibc`:
- Modifying GLT to use a 1D vector instead of 3D vector
- Test method: Modify the headers in `RangeLib.h`, run the built exectuable:
	`./range_lib -map_path ../../maps/basement_fixed.png -method glt -which_benchmark grid`
- Original timing, using `uint16_t` LUT:
	```
	...Loading range method: GiantLUTCast
	...lut size (MB): 348.129
	...construction time: 3.88341
	...Running grid benchmark
	finished grid sample after: 0.106878 sec
	-avg time per ray: 1.35079e-07 sec
	-rays cast: 676000
	-total time: 0.0913132 sec
	```
- Original timing, using `float` LUT:
    ```
	...Loading range method: GiantLUTCast
	...lut size (MB): 696.259
	...construction time: 3.84431
	...Running grid benchmark
	finished grid sample after: 0.114181 sec
	-avg time per ray: 1.4495e-07 sec
	-rays cast: 676000
	-total time: 0.0979861 sec
	```
- Modified timing, using `uint16_t` LUT:
	```
	...Loading range method: GiantLUTCast
	...lut size (MB): 348.129
	...construction time: 3.36884
	...Running grid benchmark
	finished grid sample after: 0.092248 sec
	-avg time per ray: 1.13772e-07 sec
	-rays cast: 676000
	-total time: 0.0769099 sec
	```

- Modified timing, using `float` LUT:
	```
	...Loading range method: GiantLUTCast
	...lut size (MB): 696.259
	...construction time: 3.29234
	...Running grid benchmark
	finished grid sample after: 0.112021 sec
	-avg time per ray: 1.39633e-07 sec
	-rays cast: 676000
	-total time: 0.0943918 sec
	```

As can be seen, 1D indexing is faster! It will therefore be used moving forward.

What remains to be seen is whether the `float` datatype makes a difference in actual performance.
