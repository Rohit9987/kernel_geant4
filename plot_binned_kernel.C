#include "TCanvas.h"
#include "TColor.h"
#include "TFile.h"
#include "TH2D.h"
#include "TString.h"
#include "TStyle.h"
#include "TTree.h"
#include "TGraph.h"
#include "TLine.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <iterator>
#include <limits>
#include <vector>
#include <numeric>

namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;
constexpr double kWaterDensityGramPerMm3 = 1.0e-3;
constexpr double kMeVPerGramToGy = 1.602176634e-10;

struct KernelRow {
  int radialBin{};
  int thetaBin{};
  double radialLowerMm{};
  double radialUpperMm{};
  double thetaLowerRad{};
  double thetaUpperRad{};
  double meanEdepMeVPerHistory{};
  double relativeStdError{};
};

double PositiveMinimum(const std::vector<double>& values)
{
  double minimum = std::numeric_limits<double>::infinity();
  for (const double value : values) {
    if (value > 0.0 && value < minimum) minimum = value;
  }
  return std::isfinite(minimum) ? minimum : 1.0;
}

} // namespace

void plot_binned_kernel(
  const char* inputFile = "physics/DoseKernel_1.0MeV.root",
  double displayRadiusMm = 60.0,
  const char* outputPrefix = "DoseKernel_1.0MeV")
{
  TFile file(inputFile, "READ");
  if (file.IsZombie()) {
    std::cerr << "Could not open " << inputFile << '\n';
    return;
  }

  auto* kernel = dynamic_cast<TTree*>(file.Get("Kernel"));
  if (!kernel) {
    std::cerr << "Tree 'Kernel' was not found in " << inputFile << '\n';
    file.ls();
    return;
  }

  auto* summary = dynamic_cast<TTree*>(file.Get("RunSummary"));
  if(!summary || summary->GetEntries() < 1)
  {
	  std::cerr << "Tree 'RunSummary' was not found or is empty.\n";
	  return;
  }

  if(!summary->GetBranch("Histories") ||
		  !summary->GetBranch("PrimaryEnergy_MeV") ||
		  !summary->GetBranch("LocalAtOriginEdep_MeV"))
	{
		std::cerr << "RunSummary does not contain the required branches.\n";
		summary->Print();
		return;
	}

  int histories = 0;
  double primaryEnergyMeV = 0.0;
  double localAtOriginEdepMeV = 0.0;

  summary->SetBranchAddress("Histories", &histories);
  summary->SetBranchAddress("PrimaryEnergy_MeV", & primaryEnergyMeV);
  summary->SetBranchAddress("LocalAtOriginEdep_MeV", &localAtOriginEdepMeV);

  summary->GetEntry(0);
	

  if(histories <= 0 || primaryEnergyMeV <= 0.0)
  {
	  std::cerr <<"Invalid history count or primary energy.\n";
	  return ;
  }

  const double primaryEnergyMeVPerHistory =
    primaryEnergyMeV / static_cast<double>(histories);

  const double localAtOriginMeVPerHistory =
    localAtOriginEdepMeV / static_cast<double>(histories);

  const double localAtOriginFraction =
    localAtOriginEdepMeV / primaryEnergyMeV;

  std::cout << "\nNormalization check:\n"
            << "  Primary energy per history = "
            << primaryEnergyMeVPerHistory << " MeV\n"
            << "  Local energy per history   = "
            << localAtOriginMeVPerHistory << " MeV\n"
            << "  Local energy fraction      = "
            << localAtOriginFraction << '\n';



  const char* requiredBranches[] = {
    "RadialBin",
    "ThetaBin",
    "RadialLower_mm",
    "RadialUpper_mm",
    "ThetaLower_rad",
    "ThetaUpper_rad",
    "MeanEdep_MeV_per_history",
	"RelativeStdError"

  };
  for (const char* branch : requiredBranches) {
    if (!kernel->GetBranch(branch)) {
      std::cerr << "Required branch '" << branch << "' was not found.\n";
      kernel->Print();
      return;
    }
  }

  KernelRow current;
  kernel->SetBranchAddress("RadialBin", &current.radialBin);
  kernel->SetBranchAddress("ThetaBin", &current.thetaBin);
  kernel->SetBranchAddress("RadialLower_mm", &current.radialLowerMm);
  kernel->SetBranchAddress("RadialUpper_mm", &current.radialUpperMm);
  kernel->SetBranchAddress("ThetaLower_rad", &current.thetaLowerRad);
  kernel->SetBranchAddress("ThetaUpper_rad", &current.thetaUpperRad);
  kernel->SetBranchAddress("MeanEdep_MeV_per_history",
                           &current.meanEdepMeVPerHistory);
  kernel->SetBranchAddress("RelativeStdError", &current.relativeStdError);

  std::vector<KernelRow> rows;
  rows.reserve(static_cast<std::size_t>(kernel->GetEntries()));
  int maxRadialBin = -1;
  int maxThetaBin = -1;
  for (Long64_t entry = 0; entry < kernel->GetEntries(); ++entry) {
    kernel->GetEntry(entry);
    rows.push_back(current);
    maxRadialBin = std::max(maxRadialBin, current.radialBin);
    maxThetaBin = std::max(maxThetaBin, current.thetaBin);
  }

  const int numRadialBins = maxRadialBin + 1;
  const int numThetaBins = maxThetaBin + 1;
  if (numRadialBins <= 0 || numThetaBins <= 0) {
    std::cerr << "The Kernel tree contains no valid bins.\n";
    return;
  }

  std::vector<double> radialEdges(
    static_cast<std::size_t>(numRadialBins + 1), 0.0);
  std::vector<double> thetaEdgesRad(
    static_cast<std::size_t>(numThetaBins + 1), 0.0);
  std::vector<double> doseGyPerHistory(
    static_cast<std::size_t>(numRadialBins * numThetaBins), 0.0);
  std::vector<double> relativeStdErrorPercent(
		  static_cast<std::size_t>(numRadialBins * numThetaBins), 0.0);


  // Energy deposited in each complete spherical radial shell,
  // summed over all angular sectors
  std::vector<double> radialEnergyMeVPerHistory(
		  static_cast<std::size_t>(numRadialBins), 0.0);



  for (const auto& row : rows) {
    if (row.radialBin < 0 || row.radialBin >= numRadialBins ||
        row.thetaBin < 0 || row.thetaBin >= numThetaBins) {
      continue;
    }

	radialEnergyMeVPerHistory[static_cast<std::size_t>(row.radialBin)] += 
		row.meanEdepMeVPerHistory;
    radialEdges[static_cast<std::size_t>(row.radialBin)] =
      row.radialLowerMm;
    radialEdges[static_cast<std::size_t>(row.radialBin + 1)] =
      row.radialUpperMm;
    thetaEdgesRad[static_cast<std::size_t>(row.thetaBin)] =
      row.thetaLowerRad;
    thetaEdgesRad[static_cast<std::size_t>(row.thetaBin + 1)] =
      row.thetaUpperRad;

    const double solidAngleSr = 2.0 * kPi *
      (std::cos(row.thetaLowerRad) - std::cos(row.thetaUpperRad));
    const double volumeMm3 = solidAngleSr *
      (std::pow(row.radialUpperMm, 3) -
       std::pow(row.radialLowerMm, 3)) / 3.0;
    const double massGram = volumeMm3 * kWaterDensityGramPerMm3;

    const auto linear = static_cast<std::size_t>(
      row.radialBin * numThetaBins + row.thetaBin);

	if (row.meanEdepMeVPerHistory > 0.0 &&
		std::isfinite(row.relativeStdError) &&
		row.relativeStdError >= 0.0) {

		relativeStdErrorPercent[linear] =
		100.0 * row.relativeStdError;
	}

    if (massGram > 0.0) {
      doseGyPerHistory[linear] =
        row.meanEdepMeVPerHistory / massGram * kMeVPerGramToGy;
    }
  }
  std::vector<double> cumulativeEnergyFraction(
		  static_cast<std::size_t>(numRadialBins), 0.0);
  // Begin with the energy deposited exactly at the kernel origin
  double cumulativeEnergyMeVPerHistory = localAtOriginMeVPerHistory;

  for(int radial = 0; radial < numRadialBins; ++radial)
  {
	  cumulativeEnergyMeVPerHistory += radialEnergyMeVPerHistory[static_cast<std::size_t>(radial)];
	  cumulativeEnergyFraction[static_cast<std::size_t>(radial)] = cumulativeEnergyMeVPerHistory/ primaryEnergyMeVPerHistory;
  }

    const double finalCumulativeFraction =
    cumulativeEnergyFraction.back();

  std::cout << "\nCumulative-energy check:\n"
            << "  Outer radius              = "
            << radialEdges.back() << " mm\n"
            << "  Cumulative energy fraction = "
            << finalCumulativeFraction << '\n'
            << "  Uncontained fraction       = "
            << 1.0 - finalCumulativeFraction << '\n';


  // ploting vectors
  std::vector<double> cumulativeRadiusMm(static_cast<std::size_t>(numRadialBins + 1), 0.0);
  std::vector<double> cumulativeFractionForPlot(static_cast<std::size_t>(numRadialBins + 1), 0.0);

  // The first point represents the separate deposit at r = 0;
   cumulativeRadiusMm[0] = 0.0;
   cumulativeFractionForPlot[0] = localAtOriginFraction;
  
   for (int radial = 0; radial < numRadialBins; ++radial)
   {
	   const auto plotIndex = static_cast<std::size_t>(radial+1);

	   cumulativeRadiusMm[plotIndex] = radialEdges[static_cast<std::size_t>(radial + 1)];
	   cumulativeFractionForPlot[plotIndex] = cumulativeEnergyFraction[static_cast<std::size_t>(radial)];
   }

   std::vector<double> remainingFractionForPlot(cumulativeFractionForPlot.size(), 0.0);

   for(std::size_t i = 0; i < cumulativeFractionForPlot.size(); i++)
   {
	   remainingFractionForPlot[i] = std::max(1.0e-16, 1-cumulativeFractionForPlot[i]);
   }


   struct ContainmentLevel
   {
	   const char* label;
	   double fraction;
   };

   const ContainmentLevel containmentLevels[] = {
	   {"R90", 0.90},
	   {"R95", 0.95},
	   {"R99", 0.99},
	   {"R99.9", 0.999},
	   {"R99.99", 0.9999}
   };

   std::cout << "\nEnergy-containment radii:\n";

   for(const auto& level: containmentLevels)
   {
	   const auto position = std::lower_bound(
			   cumulativeFractionForPlot.begin(),
			   cumulativeFractionForPlot.end(),
			   level.fraction);

	   if(position == cumulativeFractionForPlot.end())
	   {
		   std::cout << " " << level.label << " was not reached\n";
		   continue;
	   }

	   const std::size_t index = static_cast<std::size_t>(std::distance(
							cumulativeFractionForPlot.begin(),
							position));

	   std::cout << " " << level.label
				 << " = " << cumulativeRadiusMm[index]
				 << " mm"
				 << " (fraction = " 
				 << cumulativeFractionForPlot[index]
				 << ")\n";
   }




  const double binnedEnergyMeVPerHistory = std::accumulate(
		  radialEnergyMeVPerHistory.begin(),
		  radialEnergyMeVPerHistory.end(), 0.0);

  const double binnedEnergyFraction = 
	  binnedEnergyMeVPerHistory/ primaryEnergyMeVPerHistory;

  const double includedEnergyFraction = 
	  (binnedEnergyMeVPerHistory + localAtOriginMeVPerHistory) / primaryEnergyMeVPerHistory;

  std::cout << "\nRadial integration check:\n"
            << "  Binned energy per history = "
            << binnedEnergyMeVPerHistory << " MeV\n"
            << "  Binned energy fraction    = "
            << binnedEnergyFraction << '\n'
            << "  Including local fraction  = "
            << includedEnergyFraction << '\n'
            << "  Remaining fraction        = "
            << 1.0 - includedEnergyFraction << '\n';

  const std::vector<double> rseThresholdsPercent = {
	  1.0, 2.0, 5.0, 10.0, 20.0};

  std::vector<double> energyWithinRseThreshold(
		  rseThresholdsPercent.size(),
		  0.0);

  double energyWithInvalidRse = 0.0;
  
  for(const auto& row: rows)
  {
	  const double energy = row.meanEdepMeVPerHistory;
	  
	  if(!(energy > 0.0))
		  continue;

	  const double rsePercent = 100.0 * row.relativeStdError;

	  if(!std::isfinite(rsePercent) || rsePercent < 0.0) 
	  {
		  energyWithInvalidRse += energy;
		  continue;
	  }

	  for (std::size_t i = 0; i < rseThresholdsPercent.size(); ++i)
	  {
		  if(rsePercent <= rseThresholdsPercent[i])
		  {
			  energyWithinRseThreshold[i] += energy;
		  }
	  }
  }

  std::cout
  << "\nEnergy-weighted precision coverage:\n";

for (std::size_t i = 0;
     i < rseThresholdsPercent.size();
     ++i) {

  const double energyFraction =
    energyWithinRseThreshold[i] /
    binnedEnergyMeVPerHistory;

  std::cout
    << "  RSE <= "
    << rseThresholdsPercent[i]
    << "% : "
    << 100.0 * energyFraction
    << "% of binned energy\n";
}

std::cout
  << "  Invalid/undefined RSE : "
  << 100.0 * energyWithInvalidRse /
       binnedEnergyMeVPerHistory
  << "% of binned energy\n";


  const double maximumDose = *std::max_element(
    doseGyPerHistory.begin(), doseGyPerHistory.end());
  if (!(maximumDose > 0.0)) {
    std::cerr << "All kernel bins contain zero energy.\n";
    return;
  }

  std::vector<double> relativeDose = doseGyPerHistory;
  for (double& value : relativeDose) value /= maximumDose;

  std::vector<double> thetaEdgesDeg(thetaEdgesRad.size(), 0.0);
  for (std::size_t i = 0; i < thetaEdgesRad.size(); ++i) {
    thetaEdgesDeg[i] = thetaEdgesRad[i] * 180.0 / kPi;
  }

  gStyle->SetOptStat(0);
  gStyle->SetPalette(kViridis);

  TH2D radialAngular(
    "hKernelRTheta",
    "Volume-normalized kernel;Radius r (mm);Polar angle #theta (degrees)",
    numRadialBins, radialEdges.data(),
    numThetaBins, thetaEdgesDeg.data());
  radialAngular.SetDirectory(nullptr);

  for (int radial = 0; radial < numRadialBins; ++radial) {
    for (int theta = 0; theta < numThetaBins; ++theta) {
      const auto linear = static_cast<std::size_t>(
        radial * numThetaBins + theta);
      radialAngular.SetBinContent(radial + 1, theta + 1,
                                  relativeDose[linear]);
    }
  }

  const double relativeMinimum = std::max(
    PositiveMinimum(relativeDose), 1.0e-12);
  radialAngular.SetMinimum(relativeMinimum);
  radialAngular.SetMaximum(1.0);
  radialAngular.GetZaxis()->SetTitle("Relative dose per primary");

  TCanvas radialAngularCanvas("cKernelRTheta", "Kernel r-theta", 1000, 750);
  radialAngularCanvas.SetRightMargin(0.15);
  radialAngularCanvas.SetLogz();
  radialAngular.Draw("COLZ");
  radialAngularCanvas.SaveAs(
    TString::Format("%s_rtheta.png", outputPrefix));

  TH2D relativeUncertainty(
  "hKernelRThetaRSE",
  "Kernel statistical uncertainty;"
  "Radius r (mm);"
  "Polar angle #theta (degrees)",
  numRadialBins,
  radialEdges.data(),
  numThetaBins,
  thetaEdgesDeg.data());

	relativeUncertainty.SetDirectory(nullptr);

	for (int radial = 0; radial < numRadialBins; ++radial) {
	  for (int theta = 0; theta < numThetaBins; ++theta) {

		const auto linear = static_cast<std::size_t>(
		  radial * numThetaBins + theta);

		const double rsePercent =
		  relativeStdErrorPercent[linear];

		if (rsePercent > 0.0) {
		  relativeUncertainty.SetBinContent(
			radial + 1,
			theta + 1,
			rsePercent);
		}
	  }
	}

	const double maximumRsePercent =
  *std::max_element(
    relativeStdErrorPercent.begin(),
    relativeStdErrorPercent.end());

	const double minimumRsePercent =
	  std::max(
		PositiveMinimum(relativeStdErrorPercent),
		1.0e-3);

	relativeUncertainty.SetMinimum(minimumRsePercent);

	// Cap the displayed scale at 100% so rare, noisy bins
	// do not dominate the colour scale.
	relativeUncertainty.SetMaximum(
	  std::max(
		1.0,
		std::min(100.0, maximumRsePercent)));

	relativeUncertainty.GetZaxis()->SetTitle(
	  "Relative standard error (%)");

	TCanvas uncertaintyCanvas(
	  "cKernelRThetaRSE",
	  "Kernel statistical uncertainty",
	  1000,
	  750);

	uncertaintyCanvas.SetRightMargin(0.15);
	uncertaintyCanvas.SetLogz();

	relativeUncertainty.Draw("COLZ");

	uncertaintyCanvas.SaveAs(
	  TString::Format(
		"%s_rse_rtheta.png",
		outputPrefix));

  const int imageBins = 600;
  TH2D xz(
    "hKernelXZ",
    "Reconstructed XZ kernel;X (mm);Z (mm)",
    imageBins, -displayRadiusMm, displayRadiusMm,
    imageBins, -displayRadiusMm, displayRadiusMm);
  xz.SetDirectory(nullptr);

  std::vector<double> xzPositiveValues;
  xzPositiveValues.reserve(static_cast<std::size_t>(imageBins * imageBins));
  for (int ix = 1; ix <= imageBins; ++ix) {
    const double x = xz.GetXaxis()->GetBinCenter(ix);
    for (int iz = 1; iz <= imageBins; ++iz) {
      const double z = xz.GetYaxis()->GetBinCenter(iz);
      const double radius = std::sqrt(x * x + z * z);
      if (radius <= 0.0 || radius >= radialEdges.back()) continue;

      const auto radialUpper = std::upper_bound(
        radialEdges.begin(), radialEdges.end(), radius);
      const int radial = static_cast<int>(
        std::distance(radialEdges.begin(), radialUpper) - 1);

      const double cosTheta = std::max(-1.0, std::min(1.0, z / radius));
      const double thetaRad = std::acos(cosTheta);
      const auto thetaUpper = std::upper_bound(
        thetaEdgesRad.begin(), thetaEdgesRad.end(), thetaRad);
      int theta = static_cast<int>(
        std::distance(thetaEdgesRad.begin(), thetaUpper) - 1);
      theta = std::max(0, std::min(theta, numThetaBins - 1));

      if (radial < 0 || radial >= numRadialBins) continue;
      const auto linear = static_cast<std::size_t>(
        radial * numThetaBins + theta);
      const double value = relativeDose[linear];
      xz.SetBinContent(ix, iz, value);
      if (value > 0.0) xzPositiveValues.push_back(value);
    }
  }

  xz.SetMinimum(std::max(PositiveMinimum(xzPositiveValues), 1.0e-12));
  xz.SetMaximum(1.0);
  xz.GetZaxis()->SetTitle("Relative dose per primary");

  TCanvas xzCanvas("cKernelXZ", "Kernel XZ", 900, 850);
  xzCanvas.SetRightMargin(0.15);
  xzCanvas.SetLogz();
  xz.Draw("COLZ");
  xzCanvas.SaveAs(TString::Format("%s_xz.png", outputPrefix));


  TGraph cumulativeGraph(
		  static_cast<int>(cumulativeRadiusMm.size()),
		  cumulativeRadiusMm.data(),
		  cumulativeFractionForPlot.data());

  cumulativeGraph.SetTitle(
		  "Cumulative kernel energy;"
		  "Radius r (mm);"
		  "Cumulative fraction of primary energy");

  cumulativeGraph.SetLineColor(kBlue + 1);
  cumulativeGraph.SetLineWidth(3);
  cumulativeGraph.SetMinimum(0.0);
  cumulativeGraph.SetMaximum(1.0005);

  TCanvas cumulativeCanvas(
		  "cCumulativeEnergy",
		  "Cumulative kernel energy",
		  900,
		  700);

  cumulativeCanvas.SetGrid();
  cumulativeGraph.Draw("AL");
  cumulativeCanvas.SaveAs(
    TString::Format("%s_cumulative.png", outputPrefix));
  std::cout << "\nMaximum volume-normalized dose = " << maximumDose
            << " Gy per primary\n"
            << "Saved " << outputPrefix << "_rtheta.png, "
            << outputPrefix << "_xz.png and "
            << outputPrefix << "_cumulative.png\n";


  TGraph remainingGraph(
		  static_cast<int>(cumulativeRadiusMm.size()),
		  cumulativeRadiusMm.data(), 
		  remainingFractionForPlot.data());

  remainingGraph.SetTitle(
		  "Remaining kernel energy;"
		  "Radius r (mm);"
		  "Fraction of primary energy beyond r");

  remainingGraph.SetLineColor(kRed + 1);
  remainingGraph.SetLineWidth(3);

  const double smallestRemaining = *std::min_element(
		  remainingFractionForPlot.begin(),
		  remainingFractionForPlot.end());

  remainingGraph.SetMinimum(
		  std::max(1.0e-12, 0.5 * smallestRemaining));

  remainingGraph.SetMaximum(1.2);


  TCanvas remainingCanvas(
    "cRemainingEnergy",
    "Remaining kernel energy",
    900,
    700);

	// Keep only the vertical grid.
	// The horizontal reference lines will replace the horizontal grid.
	remainingCanvas.SetGridx();
	remainingCanvas.SetGridy(0);
	remainingCanvas.SetLogy();

   remainingGraph.Draw("AL");

  const double outerRadiusMm = cumulativeRadiusMm.back();

  TLine line99(
		  0.0, 1.0e-2, 
		  outerRadiusMm, 1.0e-2);

  TLine line999(
		  0.0, 1.0e-3,
		  outerRadiusMm, 1.0e-3);

  TLine line9999(
		  0.0, 1.0e-4,
		  outerRadiusMm, 1.0e-4);

	line99.SetLineStyle(7);
	line999.SetLineStyle(7);
	line9999.SetLineStyle(7);

	line99.SetLineWidth(3);
	line999.SetLineWidth(3);
	line9999.SetLineWidth(3);

	line99.SetLineColor(kGray + 2);
	line999.SetLineColor(kGreen + 2);
	line9999.SetLineColor(kMagenta + 2);

	line99.Draw();
	line999.Draw();
	line9999.Draw();


	remainingCanvas.Modified();
   remainingCanvas.Update();

remainingCanvas.SaveAs(
  TString::Format("%s_remaining.png", outputPrefix));


  std::cout << "\nRun summary:\n";
  summary->Show(0);

  std::cout << "\nMaximum volume-normalized dose = " << maximumDose
            << " Gy per primary\n"
            << "Saved " << outputPrefix << "_rtheta.png and "
            << outputPrefix << "_xz.png\n";
}
