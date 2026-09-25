/**
@file DistanceTransform.h
@brief Exact Euclidean distance transform of a voxel set, periodic in all
       three directions.

Answers "how deep inside the pore space does this voxel sit?", which is what
decides where chemical shrinkage empties water and where returning water
refills. A voxel at the center of a large capillary pore is many voxels from
the nearest solid; one in a narrow neck is one voxel from it. Emptying in
descending order of that depth drains the largest pores first, as the Kelvin
equation requires, and growing outward from an already-empty voxel keeps the
emptied region a compact cavity rather than a mist of isolated voxels
scattered through the paste.

Exact, not an approximation: the Felzenszwalb-Huttenlocher lower-envelope
algorithm, three O(N) passes, one per dimension. The older 7^3-window voxel
count it replaces could not distinguish a voxel deep in a large pore from one
in a small pore, because the window saturates.

Periodicity is handled by running each one-dimensional pass over the row
repeated three times and keeping the middle copy. That is exact as long as
the nearest background voxel lies within one period, which it always does:
no voxel can be more than half a box away from anything.

No THAMES dependencies. Testable in isolation; see
`src/unit_tests/test_distance_transform.cc`.
*/

#ifndef SRC_THAMESLIB_DISTANCETRANSFORM_H_
#define SRC_THAMESLIB_DISTANCETRANSFORM_H_

#include <vector>

namespace edt {

/**
@brief Get each foreground voxel's squared distance to the nearest background
voxel.

Voxel arrays are indexed `(iz * ny + iy) * nx + ix`, the order the lattice
itself uses. Distances are in voxel units, squared, so they stay exact in
floating point and can be compared without taking a root.

A background voxel gets 0. If there is no background voxel anywhere, every
voxel gets a value larger than any achievable distance, so callers that sort
by depth still behave sensibly.

@param nx is the lattice dimension along x
@param ny is the lattice dimension along y
@param nz is the lattice dimension along z
@param foreground flags the voxels to measure (size nx*ny*nz); everything
       else is background, and is what distance is measured to
@return the squared distance of every voxel, same indexing and size
*/
std::vector<double> squaredDistanceToBackground(
    int nx, int ny, int nz, const std::vector<char> &foreground);

} // namespace edt

#endif // SRC_THAMESLIB_DISTANCETRANSFORM_H_
