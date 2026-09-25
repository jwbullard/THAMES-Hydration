/**
@file Percolation.cc
@brief Method definitions for connectivity assessment.

See Percolation.h for the physics and for why each direction gets its own
labeling pass.
*/

#include "Percolation.h"

namespace percolation {

namespace {

/**
@brief Flat index of a voxel, in the same order the lattice stores sites.

@param ix is the x coordinate
@param iy is the y coordinate
@param iz is the z coordinate
@param nx is the lattice dimension along x
@param ny is the lattice dimension along y
@return the flat index (iz * ny + iy) * nx + ix
*/
inline int flatIndex(int ix, int iy, int iz, int nx, int ny) {
  return ((iz * ny) + iy) * nx + ix;
}

/**
@brief Report whether a contact between two voxels transmits load.

Both voxels are assumed to participate already. See Percolation.h for how
this reproduces VCCTL's three burn rules without naming a phase.

@param a is the flat index of the first voxel
@param b is the flat index of the second voxel
@param bonds flags voxels that grew in place (may be empty)
@param particleId holds original particle ids (may be empty)
@return true when load crosses the contact
*/
inline bool contactTransmits(int a, int b, const std::vector<char> &bonds,
                             const std::vector<int> &particleId) {
  if (!bonds.empty() && (bonds[a] || bonds[b]))
    return true;
  if (!particleId.empty() && particleId[a] > 1 &&
      particleId[a] == particleId[b])
    return true;
  return false;
}

/**
@brief Assess one axis with that axis non-periodic and the other two periodic.

Labels every participating voxel by flood fill, recording for each cluster
whether it touches the low face and the high face of the spanning axis. A
cluster touching both is a through path.

@param axis is 0 for x, 1 for y, 2 for z
@param nx is the lattice dimension along x
@param ny is the lattice dimension along y
@param nz is the lattice dimension along z
@param participates flags each voxel as part of the burning set
@param bonds flags each voxel as bonding to whatever it touches
@param particleId holds original particle ids
@param numParticipating is the total size of the burning set
@param label is scratch space of size nx*ny*nz, overwritten here
@param queue is scratch space for the flood fill, overwritten here
@return the spanning flag, spanning-voxel count, and connected fraction
*/
DirectionResult assessAxis(int axis, int nx, int ny, int nz,
                           const std::vector<char> &participates,
                           const std::vector<char> &bonds,
                           const std::vector<int> &particleId,
                           int numParticipating, std::vector<int> &label,
                           std::vector<int> &queue) {
  DirectionResult result;

  const int dim[3] = {nx, ny, nz};
  const int numSites = nx * ny * nz;

  label.assign(numSites, -1);

  int nextLabel = 0;
  int numSpanning = 0;

  for (int seed = 0; seed < numSites; ++seed) {
    if (!participates[seed] || label[seed] >= 0)
      continue;

    // Flood fill this cluster, noting which faces of the spanning axis it
    // reaches and how many voxels it holds.
    const int thisLabel = nextLabel++;
    bool touchesLow = false;
    bool touchesHigh = false;
    int clusterSize = 0;

    queue.clear();
    queue.push_back(seed);
    label[seed] = thisLabel;

    for (std::size_t head = 0; head < queue.size(); ++head) {
      const int current = queue[head];
      ++clusterSize;

      const int ix = current % nx;
      const int iy = (current / nx) % ny;
      const int iz = current / (nx * ny);
      const int coord[3] = {ix, iy, iz};

      if (coord[axis] == 0)
        touchesLow = true;
      if (coord[axis] == dim[axis] - 1)
        touchesHigh = true;

      // Six face neighbours. The spanning axis does not wrap, so a path
      // cannot percolate by leaving one face and re-entering the opposite
      // one; the other two axes do wrap, because the microstructure is a
      // periodic RVE and a load path may legitimately cross the seam there.
      for (int d = 0; d < 3; ++d) {
        for (int step = -1; step <= 1; step += 2) {
          int neighbor[3] = {ix, iy, iz};
          neighbor[d] = coord[d] + step;

          if (neighbor[d] < 0) {
            if (d == axis)
              continue;
            neighbor[d] += dim[d];
          } else if (neighbor[d] >= dim[d]) {
            if (d == axis)
              continue;
            neighbor[d] -= dim[d];
          }

          const int neighborIndex =
              flatIndex(neighbor[0], neighbor[1], neighbor[2], nx, ny);
          if (!participates[neighborIndex] || label[neighborIndex] >= 0)
            continue;
          if (!contactTransmits(current, neighborIndex, bonds, particleId))
            continue;

          label[neighborIndex] = thisLabel;
          queue.push_back(neighborIndex);
        }
      }
    }

    if (touchesLow && touchesHigh) {
      result.spans = true;
      numSpanning += clusterSize;
    }
  }

  result.numSpanningVoxels = numSpanning;
  result.connectedFraction =
      (numParticipating > 0)
          ? static_cast<double>(numSpanning) / static_cast<double>(numParticipating)
          : 0.0;
  return result;
}

} // namespace

Result assess(int nx, int ny, int nz, const std::vector<char> &participates,
              const std::vector<char> &bonds,
              const std::vector<int> &particleId) {
  Result result;

  const int numSites = nx * ny * nz;
  if (numSites <= 0 || static_cast<int>(participates.size()) != numSites)
    return result;

  for (int i = 0; i < numSites; ++i) {
    if (participates[i])
      ++result.numParticipating;
  }
  if (result.numParticipating == 0)
    return result;

  // Scratch space reused across the three passes.
  std::vector<int> label;
  std::vector<int> queue;
  queue.reserve(result.numParticipating);

  for (int axis = 0; axis < 3; ++axis) {
    result.direction[axis] =
        assessAxis(axis, nx, ny, nz, participates, bonds, particleId,
                   result.numParticipating, label, queue);
  }

  return result;
}

} // namespace percolation
