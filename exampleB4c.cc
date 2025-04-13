#include "DetectorConstruction.hh"
#include "ActionInitialization.hh"
#include "PhysicsList.hh"

#include "G4RunManagerFactory.hh"
#include "G4SteppingVerbose.hh"
#include "G4UIcommand.hh"
#include "G4UImanager.hh"
#include "G4UIExecutive.hh"
#include "G4VisExecutive.hh"
#include "Randomize.hh"

#include "FTFP_BERT.hh"


// Usage message
namespace {
  void PrintUsage() {
    G4cerr << " Usage: " << G4endl;
    G4cerr << " exampleB4c [-m macro] [-e energy] [-u UIsession] [-t nThreads] [-vDefault]" << G4endl;
    G4cerr << "   -e energy : photon energy (e.g., 0.1MeV, 1.0MeV)" << G4endl;
    G4cerr << "   note: -t option is available only for multi-threaded mode." << G4endl;
  }
}

int main(int argc,char** argv)
{
  if ( argc > 9 ) {
    PrintUsage();
    return 1;
  }
  // random

  G4String macro;
  G4String session;
  G4String energyStr = "1.0MeV"; // Default photon energy
  G4bool verboseBestUnits = true;

#ifdef G4MULTITHREADED
  G4int nThreads = 0;
#endif

  // Parse command-line arguments
  for (G4int i=1; i<argc; i+=2) {
    if      ( G4String(argv[i]) == "-m" ) macro = argv[i+1];
    else if ( G4String(argv[i]) == "-u" ) session = argv[i+1];
    else if ( G4String(argv[i]) == "-e" ) energyStr = argv[i+1];
#ifdef G4MULTITHREADED
    else if ( G4String(argv[i]) == "-t" ) nThreads = G4UIcommand::ConvertToInt(argv[i+1]);
#endif
    else if ( G4String(argv[i]) == "-vDefault" ) {
      verboseBestUnits = false;
      --i;  // no parameter follows this option
    }
    else {
      PrintUsage();
      return 1;
    }
  }

  // Interactive mode detection
  G4UIExecutive* ui = nullptr;
  if ( ! macro.size() ) {
    ui = new G4UIExecutive(argc, argv, session);
  }

  if ( verboseBestUnits ) {
    G4SteppingVerbose::UseBestUnit(4);
  }

  // Construct run manager
  auto* runManager =
    G4RunManagerFactory::CreateRunManager(G4RunManagerType::Default);

#ifdef G4MULTITHREADED
  if ( nThreads > 0 ) runManager->SetNumberOfThreads(nThreads);
#endif

  // Mandatory initialization
  runManager->SetUserInitialization(new B4c::DetectorConstruction());
  runManager->SetUserInitialization(new B4::PhysicsList);

  // Pass energy string to ActionInitialization
  runManager->SetUserInitialization(new B4c::ActionInitialization(energyStr));

  // Visualization manager
  auto visManager = new G4VisExecutive;
  visManager->Initialize();

  // UI Manager
  auto UImanager = G4UImanager::GetUIpointer();

  // Execute macro if provided
  if ( macro.size() ) {
    UImanager->ApplyCommand("/control/execute " + macro);
  } else {
    UImanager->ApplyCommand("/control/execute init_vis.mac");
    if (ui->IsGUI()) UImanager->ApplyCommand("/control/execute gui.mac");
    ui->SessionStart();
    delete ui;
  }

  delete visManager;
  delete runManager;
}

