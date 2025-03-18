#include "DetectorConstruction.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"

#include "G4SDManager.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4LogicalVolumeStore.hh"
#include "VoxelSensitiveDetector.hh"

#include "G4UserLimits.hh"


namespace B4c
{


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreadLocal
G4GlobalMagFieldMessenger* DetectorConstruction::fMagFieldMessenger = nullptr;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
DetectorConstruction::DetectorConstruction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  // Define materials
  DefineMaterials();

  // Define volumes
  return DefineVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{
  // Lead material defined using NIST Manager
  auto nistManager = G4NistManager::Instance();

  nistManager->FindOrBuildMaterial("G4_WATER");
  // Print materials
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::DefineVolumes()
{
  G4NistManager* nist = G4NistManager::Instance();
  G4Material* water = nist->FindOrBuildMaterial("G4_WATER");

	G4double worldSizeXY = 0.5*m,
			 worldSizeZ = 0.5*m;
  //
  // World
  //
  auto worldS
    = new G4Box("World",           // its name
                 worldSizeXY/2, worldSizeXY/2, worldSizeZ/2); // its size

  auto worldLV
    = new G4LogicalVolume(
                 worldS,           // its solid
                 water,  // its material
                 "World");         // its name

  auto worldPV
    = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(),  // at (0,0,0)
                 worldLV,          // its logical volume
                 "World",          // its name
                 0,                // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps


	// voxels
	// definitions
	G4double halfVoxelRes = 0.5*mm;
	G4int numVoxels = 101;

	G4Box* voxelBox = new G4Box("Voxel", halfVoxelRes, halfVoxelRes, halfVoxelRes);
    G4LogicalVolume* voxelLV = new G4LogicalVolume(voxelBox, water, "VoxelLV");

	G4double startX = -(numVoxels - 1) * halfVoxelRes;
	G4double startY = -(numVoxels - 1) * halfVoxelRes;
	G4double startZ = -(numVoxels - 1) * halfVoxelRes;

	for (int i = 0; i < numVoxels; ++i) {
		for (int j = 0; j < numVoxels; ++j) {
			for (int k = 0; k < numVoxels; ++k) {

				G4double xPos = startX + i*2*halfVoxelRes;
				G4double yPos = startY + j*2*halfVoxelRes;
				G4double zPos = startZ + k*2*halfVoxelRes;

				G4ThreeVector position(xPos, yPos, zPos);
				G4int id = i*numVoxels*numVoxels + j*numVoxels + k;

				new G4PVPlacement(nullptr, position, voxelLV, "Voxel", worldLV, false, id);
			}
		}
	} 

	// Define a step limit (in mm)
    G4double maxStep = 0.1 * mm;  
    G4UserLimits* stepLimit = new G4UserLimits(maxStep);

    // Get the logical volume you want to apply the limit to
    voxelLV->SetUserLimits(stepLimit);



  //
  // Always return the physical World
  //
  return worldPV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructSDandField()
{
        G4SDManager* sdManager = G4SDManager::GetSDMpointer();
        auto* voxelSD = new VoxelSensitiveDetector("VoxelSD");
        G4LogicalVolume* voxelLV = G4LogicalVolumeStore::GetInstance()->GetVolume("VoxelLV");
        if (voxelLV) voxelLV->SetSensitiveDetector(voxelSD);
        sdManager->AddNewDetector(voxelSD);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}
