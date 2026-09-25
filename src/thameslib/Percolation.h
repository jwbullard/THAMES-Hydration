/**
@file Percolation.h
@brief Connectivity assessment of a voxel set, for set detection and for
       capillary-pore percolation.

Port of VCCTL's `burnset` and `burn3d` (disrealnew), generalized so that no
phase name appears anywhere: which voxels take part, and which contacts
transmit load, arrive as plain per-voxel arrays. The same function therefore
serves both uses.

**Rigidity percolation (set detection).** `participates` marks solid voxels
whose phase declares `rigidity.participates` in simparams.json; `bonds` marks
voxels that grew or transformed in place, which adhere to whatever they
touch; `particleId` carries the original particle ids read from the `.pimg`
image. A contact transmits load when

    bonds[a] || bonds[b] || (particleId[a] > 1 && particleId[a] == particleId[b])

reproducing VCCTL's three burn rules: burning spreads into a bridging
hydrate, from a hydrate into a reactant, and between two reactant voxels only
inside one original grain. That last clause is what stops grains which merely
touch after flocculation from reading as a connected solid path, which would
otherwise report set at t = 0.

**Capillary percolation.** `participates` marks VOID and ELECTROLYTE voxels,
`bonds` is all true, and `particleId` is empty: connectivity is then purely
geometric.

**Periodicity.** Each direction is assessed in its own labeling pass with
that axis NON-periodic and the other two periodic, as VCCTL did. A path may
legitimately wrap in y and z while spanning x — the microstructure is a
periodic RVE — but a path wrapping along the spanning axis itself would
percolate trivially. A single fully-periodic labeling followed by a spanning
test is exactly the bug fixed in the UI connectivity calculator in September
2026: the periodic merge stitched isolated grain cores across the seam and
then declared percolation on the merged label.

Costs three passes over the lattice, tens of milliseconds at 100^3, and runs
only on scheduled cycles (every 10 min of hydration time before set, every
1-2 h after).

No THAMES dependencies. Testable in isolation; see
`src/unit_tests/test_percolation.cc`.
*/

#ifndef SRC_THAMESLIB_PERCOLATION_H_
#define SRC_THAMESLIB_PERCOLATION_H_

#include <vector>

namespace percolation {

/**
@struct DirectionResult
@brief Outcome of assessing one axis.
*/
struct DirectionResult {
  bool spans = false; /**< a cluster reaches both faces along this axis */
  int numSpanningVoxels = 0;       /**< participating voxels in clusters that
                                        touch both faces */
  double connectedFraction = 0.0;  /**< numSpanningVoxels / numParticipating,
                                        0 when nothing participates */
};

/**
@struct Result
@brief Outcome of a complete assessment, one entry per axis.
*/
struct Result {
  DirectionResult direction[3]; /**< indexed 0 = x, 1 = y, 2 = z */
  int numParticipating = 0;     /**< voxels in the burning set */

  /**
  @brief Report whether the set spans all three directions.

  This is the initial-set criterion: rigidity percolation in 3D, the first
  moment a load path exists along every axis.

  @return true when every direction spans
  */
  bool spansAllDirections() const {
    return (direction[0].spans && direction[1].spans && direction[2].spans);
  }

  /**
  @brief Get the smallest connected fraction over the three directions.

  Used for the final-set criterion, which compares against a threshold held
  by the caller (0.985, VCCTL's value). Taking the minimum keeps final set
  consistent with initial set in demanding all three directions rather than
  just the most favourable one.

  @return the smallest per-direction connected fraction
  */
  double minConnectedFraction() const {
    double smallest = direction[0].connectedFraction;
    for (int d = 1; d < 3; ++d) {
      if (direction[d].connectedFraction < smallest)
        smallest = direction[d].connectedFraction;
    }
    return smallest;
  }
};

/**
@brief Assess how well the participating voxel set is connected.

Voxel arrays are indexed `(iz * ny + iy) * nx + ix`, the order the lattice
itself uses. `bonds` and `particleId` may be empty: an empty `bonds` means no
voxel bonds on contact by virtue of having grown in place, and an empty
`particleId` disables the same-grain rule. Passing `bonds` all true with an
empty `particleId` gives plain geometric connectivity, which is what
capillary percolation wants.

@param nx is the lattice dimension along x
@param ny is the lattice dimension along y
@param nz is the lattice dimension along z
@param participates flags each voxel as part of the burning set (size nx*ny*nz)
@param bonds flags each voxel as bonding to whatever it touches (size nx*ny*nz
       or empty)
@param particleId is the original particle id of each voxel (size nx*ny*nz or
       empty); values <= 1 mean "not part of any original particle"
@return the per-direction spanning flags and connected fractions
*/
Result assess(int nx, int ny, int nz, const std::vector<char> &participates,
              const std::vector<char> &bonds,
              const std::vector<int> &particleId);

} // namespace percolation

#endif // SRC_THAMESLIB_PERCOLATION_H_
