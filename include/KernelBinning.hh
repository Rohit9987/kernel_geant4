#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

namespace B4c {
namespace KernelBinning {

constexpr double kPi = 3.141592653589793238462643383279502884;
constexpr std::size_t kNumThetaBins = 48;
constexpr int kInvalidBin = -1;

// Radial resolution is deliberately finest over the charged-particle core
// and progressively coarser in the scattered-photon tail.
inline const std::vector<double>& RadialEdgesMm()
{
  static const std::vector<double> edges = [] {
    std::vector<double> values;

    const auto appendUniformRange = [&values](double lower,
                                               double upper,
                                               double width) {
      if (values.empty()) values.push_back(lower);
      const auto nBins = static_cast<std::size_t>(
        std::llround((upper - lower) / width));
      for (std::size_t i = 1; i <= nBins; ++i) {
        values.push_back(lower + static_cast<double>(i) * width);
      }
    };

    appendUniformRange(0.0,   10.0,  0.1);
    appendUniformRange(10.0,  50.0,  0.5);
    appendUniformRange(50.0, 200.0,  2.0);
    appendUniformRange(200.0, 600.0, 10.0);
    appendUniformRange(600.0, 1000.0, 20.0);
    appendUniformRange(1000.0, 1500.0, 50.0);
    appendUniformRange(1500.0, 2500.0, 100.0);
    appendUniformRange(2500.0, 3000.0, 200.0);

    return values;
  }();
  return edges;
}

inline std::size_t NumRadialBins()
{
  return RadialEdgesMm().size() - 1;
}

inline std::size_t NumBins()
{
  return NumRadialBins() * kNumThetaBins;
}

inline double ThetaWidthRad()
{
  return kPi / static_cast<double>(kNumThetaBins);
}

inline int FindRadialBin(double radiusMm)
{
  const auto& edges = RadialEdgesMm();
  if (!std::isfinite(radiusMm) || radiusMm < edges.front() ||
      radiusMm >= edges.back()) {
    return kInvalidBin;
  }

  const auto upper = std::upper_bound(edges.begin(), edges.end(), radiusMm);
  return static_cast<int>(std::distance(edges.begin(), upper) - 1);
}

inline int FindThetaBin(double thetaRad)
{
  if (!std::isfinite(thetaRad) || thetaRad < 0.0 || thetaRad > kPi) {
    return kInvalidBin;
  }
  if (thetaRad == kPi) return static_cast<int>(kNumThetaBins - 1);

  const auto index = static_cast<std::size_t>(thetaRad / ThetaWidthRad());
  return static_cast<int>(std::min(index, kNumThetaBins - 1));
}

inline std::size_t LinearIndex(std::size_t radialIndex,
                               std::size_t thetaIndex)
{
  return radialIndex * kNumThetaBins + thetaIndex;
}

inline std::size_t RadialIndex(std::size_t linearIndex)
{
  return linearIndex / kNumThetaBins;
}

inline std::size_t ThetaIndex(std::size_t linearIndex)
{
  return linearIndex % kNumThetaBins;
}

inline double RadialLowerMm(std::size_t radialIndex)
{
  return RadialEdgesMm().at(radialIndex);
}

inline double RadialUpperMm(std::size_t radialIndex)
{
  return RadialEdgesMm().at(radialIndex + 1);
}

inline double ThetaLowerRad(std::size_t thetaIndex)
{
  return static_cast<double>(thetaIndex) * ThetaWidthRad();
}

inline double ThetaUpperRad(std::size_t thetaIndex)
{
  return static_cast<double>(thetaIndex + 1) * ThetaWidthRad();
}

inline double SolidAngleSr(std::size_t thetaIndex)
{
  return 2.0 * kPi *
    (std::cos(ThetaLowerRad(thetaIndex)) -
     std::cos(ThetaUpperRad(thetaIndex)));
}

inline double BinVolumeMm3(std::size_t radialIndex,
                           std::size_t thetaIndex)
{
  const double r0 = RadialLowerMm(radialIndex);
  const double r1 = RadialUpperMm(radialIndex);
  return SolidAngleSr(thetaIndex) * (r1 * r1 * r1 - r0 * r0 * r0) / 3.0;
}

} // namespace KernelBinning
} // namespace B4c
