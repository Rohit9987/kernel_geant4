#pragma once

#include <cstddef>

namespace B4c {

// Energy deposits that cannot be assigned to a radial-angular kernel bin.
// Keeping the reasons separate makes energy-closure and geometry checks
// possible without writing individual Geant4 steps to the ROOT file.
enum class UnbinnedReason : std::size_t {
  LocalAtOrigin = 0,
  OutsideKernelRadius,
  MissingKernelFrame,
  InvalidDirectionOrAngle,
  Count
};

constexpr std::size_t kNumUnbinnedReasons =
  static_cast<std::size_t>(UnbinnedReason::Count);

constexpr std::size_t UnbinnedReasonIndex(UnbinnedReason reason)
{
  return static_cast<std::size_t>(reason);
}

static_assert(kNumUnbinnedReasons == 4,
              "Update all diagnostic energy counters when adding a reason.");

} // namespace B4c
