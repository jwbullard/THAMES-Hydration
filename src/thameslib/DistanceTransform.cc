/**
@file DistanceTransform.cc
@brief Method definitions for the exact Euclidean distance transform.

See DistanceTransform.h for what the transform is used for and why it is
exact.
*/

#include "DistanceTransform.h"

#include <vector>

namespace edt {

namespace {

/**
@brief Transform one row by the lower envelope of its parabolas.

The Felzenszwalb-Huttenlocher step. Each sample f[q] is the vertex of a
parabola (x - q)^2 + f[q]; the transform of the row is the lower envelope of
all of them, evaluated at every sample. Walking left to right and keeping
only the parabolas still visible makes this linear in the row length.

@param f is the row to transform
@param d receives the transformed row, resized as needed
@param v is scratch space for the visible parabolas' vertices
@param z is scratch space for the boundaries between them
*/
void lowerEnvelope(const std::vector<double> &f, std::vector<double> &d,
                   std::vector<int> &v, std::vector<double> &z) {
  const int n = static_cast<int>(f.size());
  d.resize(n);
  v.resize(n);
  z.resize(n + 1);

  const double kFar = 1.0e300;

  int k = 0;
  v[0] = 0;
  z[0] = -kFar;
  z[1] = kFar;

  for (int q = 1; q < n; ++q) {
    // Where does q's parabola cross the one currently on top? If that is to
    // the left of where the top one became visible, the top one is hidden
    // from here on and must be dropped.
    double s = 0.0;
    while (true) {
      const int p = v[k];
      s = ((f[q] + static_cast<double>(q) * q) -
           (f[p] + static_cast<double>(p) * p)) /
          (2.0 * static_cast<double>(q - p));
      if (s > z[k])
        break;
      if (k == 0) {
        s = -kFar;
        break;
      }
      --k;
    }
    ++k;
    v[k] = q;
    z[k] = s;
    z[k + 1] = kFar;
  }

  k = 0;
  for (int q = 0; q < n; ++q) {
    while (z[k + 1] < static_cast<double>(q))
      ++k;
    const double dx = static_cast<double>(q - v[k]);
    d[q] = dx * dx + f[v[k]];
  }
}

/**
@brief Transform one row as if it repeated forever.

Runs the envelope over three copies of the row and keeps the middle one, so a
voxel near one end sees the far end as its neighbour. Exact because nothing
can be more than half a period from its nearest background voxel, which is
well inside the copies on either side.

@param row is the row to transform, overwritten with the result
@param tiled is scratch space for the three copies
@param d is scratch space for the transformed copies
@param v is scratch space for the visible parabolas' vertices
@param z is scratch space for the boundaries between them
*/
void lowerEnvelopePeriodic(std::vector<double> &row,
                           std::vector<double> &tiled, std::vector<double> &d,
                           std::vector<int> &v, std::vector<double> &z) {
  const int n = static_cast<int>(row.size());
  if (n <= 0)
    return;

  tiled.resize(3 * n);
  for (int i = 0; i < n; ++i) {
    tiled[i] = tiled[i + n] = tiled[i + 2 * n] = row[i];
  }

  lowerEnvelope(tiled, d, v, z);

  for (int i = 0; i < n; ++i)
    row[i] = d[i + n];
}

} // namespace

std::vector<double> squaredDistanceToBackground(
    int nx, int ny, int nz, const std::vector<char> &foreground) {
  std::vector<double> dist;

  const int numSites = nx * ny * nz;
  if (numSites <= 0 || static_cast<int>(foreground.size()) != numSites)
    return dist;

  // Larger than any squared distance the box can hold, but small enough to
  // add to without trouble. Used both to seed the foreground and as the
  // answer when there is no background at all.
  const double kUnreachable =
      4.0 * (static_cast<double>(nx) * nx + static_cast<double>(ny) * ny +
             static_cast<double>(nz) * nz);

  dist.resize(numSites);
  for (int i = 0; i < numSites; ++i)
    dist[i] = foreground[i] ? kUnreachable : 0.0;

  // Scratch reused by every row of every pass.
  std::vector<double> row, tiled, d, z;
  std::vector<int> v;

  // Pass along x: rows are contiguous.
  row.resize(nx);
  for (int iz = 0; iz < nz; ++iz) {
    for (int iy = 0; iy < ny; ++iy) {
      const int base = (iz * ny + iy) * nx;
      for (int ix = 0; ix < nx; ++ix)
        row[ix] = dist[base + ix];
      lowerEnvelopePeriodic(row, tiled, d, v, z);
      for (int ix = 0; ix < nx; ++ix)
        dist[base + ix] = row[ix];
    }
  }

  // Pass along y.
  row.resize(ny);
  for (int iz = 0; iz < nz; ++iz) {
    for (int ix = 0; ix < nx; ++ix) {
      for (int iy = 0; iy < ny; ++iy)
        row[iy] = dist[(iz * ny + iy) * nx + ix];
      lowerEnvelopePeriodic(row, tiled, d, v, z);
      for (int iy = 0; iy < ny; ++iy)
        dist[(iz * ny + iy) * nx + ix] = row[iy];
    }
  }

  // Pass along z.
  row.resize(nz);
  for (int iy = 0; iy < ny; ++iy) {
    for (int ix = 0; ix < nx; ++ix) {
      for (int iz = 0; iz < nz; ++iz)
        row[iz] = dist[(iz * ny + iy) * nx + ix];
      lowerEnvelopePeriodic(row, tiled, d, v, z);
      for (int iz = 0; iz < nz; ++iz)
        dist[(iz * ny + iy) * nx + ix] = row[iz];
    }
  }

  return dist;
}

} // namespace edt
