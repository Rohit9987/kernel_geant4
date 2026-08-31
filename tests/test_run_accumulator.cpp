#include "Run.hh"

#include <array>
#include <cmath>
#include <cstddef>
#include <iostream>

namespace {

int failures = 0;

void Check(bool condition, const char* expression, int line)
{
  if (condition) return;

  std::cerr << "FAILED at line " << line << ": " << expression << '\n';
  ++failures;
}

bool NearlyEqual(double a, double b, double tolerance = 1.0e-12)
{
  return std::abs(a - b) <= tolerance;
}

} // namespace

#define CHECK(expression) Check((expression), #expression, __LINE__)

int main()
{
  MyRun run;
  run.AddPhoton();
  run.AddPhoton();

  run.AddEventBin(0, 0.25);
  run.AddEventBin(0, 0.75);

  std::array<G4double, B4c::kNumUnbinnedReasons> eventOne{};
  eventOne[B4c::UnbinnedReasonIndex(
    B4c::UnbinnedReason::LocalAtOrigin)] = 1.0;
  eventOne[B4c::UnbinnedReasonIndex(
    B4c::UnbinnedReason::OutsideKernelRadius)] = 2.0;
  run.AddUnbinnedEvent(eventOne);

  std::array<G4double, B4c::kNumUnbinnedReasons> eventTwo{};
  eventTwo[B4c::UnbinnedReasonIndex(
    B4c::UnbinnedReason::LocalAtOrigin)] = 3.0;
  eventTwo[B4c::UnbinnedReasonIndex(
    B4c::UnbinnedReason::MissingKernelFrame)] = 4.0;
  eventTwo[B4c::UnbinnedReasonIndex(
    B4c::UnbinnedReason::InvalidDirectionOrAngle)] = 5.0;
  run.AddUnbinnedEvent(eventTwo);

  run.AddEscapedEvent(0.5);
  run.AddEscapedEvent(1.5);
  run.AddPrimaryEnergyEvent(1.0);
  run.AddPrimaryEnergyEvent(1.0);

  CHECK(run.GetPhotonCount() == 2);
  CHECK(NearlyEqual(run.GetEdepSum(0), 1.0));
  CHECK(NearlyEqual(run.GetEdepSumSquares(0), 0.625));
  CHECK(NearlyEqual(run.GetTotalScoredEdep(), 1.0));

  CHECK(NearlyEqual(run.GetUnbinnedEdepSum(
    B4c::UnbinnedReason::LocalAtOrigin), 4.0));
  CHECK(NearlyEqual(run.GetUnbinnedEdepSumSquares(
    B4c::UnbinnedReason::LocalAtOrigin), 10.0));
  CHECK(NearlyEqual(run.GetUnbinnedEdepSum(
    B4c::UnbinnedReason::OutsideKernelRadius), 2.0));
  CHECK(NearlyEqual(run.GetUnbinnedEdepSum(
    B4c::UnbinnedReason::MissingKernelFrame), 4.0));
  CHECK(NearlyEqual(run.GetUnbinnedEdepSum(
    B4c::UnbinnedReason::InvalidDirectionOrAngle), 5.0));
  CHECK(NearlyEqual(run.GetTotalUnbinnedEdepSum(), 15.0));
  CHECK(NearlyEqual(run.GetTotalUnbinnedEdepSumSquares(), 153.0));

  CHECK(NearlyEqual(run.GetEscapedEnergySum(), 2.0));
  CHECK(NearlyEqual(run.GetEscapedEnergySumSquares(), 2.5));
  CHECK(NearlyEqual(run.GetPrimaryEnergySum(), 2.0));
  CHECK(NearlyEqual(run.GetPrimaryEnergySumSquares(), 2.0));

  MyRun worker;
  worker.AddPhoton();
  worker.AddEventBin(0, 2.0);
  worker.AddEscapedEvent(0.25);
  worker.AddPrimaryEnergyEvent(1.0);

  std::array<G4double, B4c::kNumUnbinnedReasons> workerEvent{};
  workerEvent[B4c::UnbinnedReasonIndex(
    B4c::UnbinnedReason::OutsideKernelRadius)] = 0.5;
  worker.AddUnbinnedEvent(workerEvent);

  run.Merge(&worker);

  CHECK(run.GetPhotonCount() == 3);
  CHECK(NearlyEqual(run.GetEdepSum(0), 3.0));
  CHECK(NearlyEqual(run.GetEdepSumSquares(0), 4.625));
  CHECK(NearlyEqual(run.GetTotalScoredEdep(), 3.0));
  CHECK(NearlyEqual(run.GetUnbinnedEdepSum(
    B4c::UnbinnedReason::OutsideKernelRadius), 2.5));
  CHECK(NearlyEqual(run.GetTotalUnbinnedEdepSum(), 15.5));
  CHECK(NearlyEqual(run.GetTotalUnbinnedEdepSumSquares(), 153.25));
  CHECK(NearlyEqual(run.GetEscapedEnergySum(), 2.25));
  CHECK(NearlyEqual(run.GetPrimaryEnergySum(), 3.0));

  if (failures != 0) {
    std::cerr << failures << " run-accumulator test(s) failed.\n";
    return 1;
  }

  std::cout << "Run accumulator tests passed\n";
  return 0;
}

#undef CHECK
