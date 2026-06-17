"""
This script is used to run the GATE code using multiple cores.  It creates macro scirpts where time has
been calcualted for each core.  Since each core processes a different macro file they will each have a
separate *.root file.  These separate root files are combined together using hadd which combines the 
leaves of the multiple output *.root files into a single leaf which is then stored in a root file.  

Author: Rayhaan Perin
Last Revised: 15-09-2021
Required Packages: numpy

To install packages packages use pip3 install package or pip install package,
alternatively if one uses conda you can run conda install package.  Conda blows up G4.  
"""

######################
###### PACKAGES ######
######################

import os 
import multiprocessing
import numpy as np
import shutil
from subprocess import run, PIPE
import re

###########################
###### FILE CONTROLS ######
###########################

def createFolderIfNotExist(folder):
    if not os.path.exists(folder):
            os.makedirs(folder)

radii = [5, 4, 3, 2, 1, 0.5, 0.1]
tangentialVelocity = 10 # mm/s

for h in range(len(radii)):
    
    period = (2*np.pi*radii[h])/tangentialVelocity
    totalTime = 2*period

    NAME = "GATE_Multiprocess" # Doesn't do anything anymore 
    DIR = "/home/rayhaan/REPO_HR++/GATE_HR/" # Doesn't do anything anymore
    MACRO_FOLDER = "/home/rayhaan/REPO_HR++/GATE_HR/MacroDir/" # The folder where to store the macro files for multiprocessing, temporary, will be removed afterwards
    ROOT_FOLDER = "/home/rayhaan/REPO_HR++/GATE_HR/processOut/" # The folder to store the output root files from multiprocessing, temporary, will be removed afterwards
    OUTPUT_ROOT = "/home/rayhaan/REPO_HR++/GATE_HR/Circles/ROOT/" # Folder to write final output root file
    OUTPUT_ROOT_FILE_NAME = "Circle_radius_{}mm.root".format(radii[h]) # The output root file name 

    createFolderIfNotExist(MACRO_FOLDER)
    createFolderIfNotExist(ROOT_FOLDER)

    ###########################
    ###### CORE CONTROLS ######
    ###########################

    cores = 50 # How many processes to create, speeds up sim 

    ###########################
    ###### MACRO CONTROLS #####
    ###########################

    # TRACER SIZE
    tracerSize = 150 # Tracer Radius 

    #   WORLD
    worldXLength = "1.1 m" 
    worldYLength = "1.1 m"
    worldZLength = "30. cm"
    worldMaterial = "Air"

    #   PEPT GEOMETRY
    # CPET
    CPETcylinderMaterial = "Air"
    CPETcylinderMaxRadius = "500.0 mm"
    CPETcylinderMinRadius = "420.0 mm"
    CPETclyinderHeight = "23.4 cm"

    # Sector
    sectorMaterial = "Air"
    sectorTranslation = "0.0 0.0 0.0 mm"
    sectorMaxRadius = "500.0 mm"
    sectorMinRadius = "420.0 mm"
    sectorHeight = "23.4 cm"
    sectorPhiStart = "-5.0 deg"
    sectorDeltaPhi = "5 deg"

    # Sector Repeater
    sectorRepeatType = "ring"
    sectorRepeatNo = "72"
    # sectorRepeatExclude1 = "7"   # To model the missing buckets in the validation data.  Leave in for general use.  
    # sectorRepeatExclude2 = "11"  #

    # Casette
    cassetteMaterial = "Air"
    cassetteTranslation = "0.0 0.0 0.0 mm"
    cassetteMaxRadius = "500.0 mm"
    cassetteMinRadius = "420.0 mm"
    cassetteHeight = "23.4 cm"
    cassettePhiStart = "-5.0 deg"
    cassetteDeltaPhi = "5 deg"

    # block 
    blockMaterial = "Air"
    blockTranslation = "459.5621819276546 -20.06491818805456 0.0 mm"
    blockRotationAxis = "0 0 1"
    blockRotationAngle = "-2.5 deg"
    blockXLength = "80.0 mm"
    blockYLength = "35.62 mm"
    blockZLength = "38.34 mm"

    # block repeat
    blockRepeatType = "linear"
    blockRepeatNo = "6"
    blockRepeatVector = "0.0 0.0 39.132 mm"

    # crystal 
    crystalMaterial = "BGO"
    crystalXLength = "30.0 mm"
    crystalYLength = "4.05 mm"
    crystalZLength = "4.39 mm"
    crystalTranslation = "-25.0 0.0 0.0 mm"

    # crystal repeat
    crystalRepeatType = "cubicArray"
    crystalRepeatX = "1"
    crystalRepeatY = "8"
    crystalRepeatZ = crystalRepeatY
    crystalRepeatVector = "0. 4.61 4.85 mm"

    # PMT outer layer
    pmtOuterMaterial = "Borosilicate-Glass"
    pmtOuterTranslation = "15.0 0.0 0.0 mm"
    pmtOuterRadiusMax = "8.79 mm"
    pmtOuterRadiusMin = "7.79 mm"
    pmtOuterHeight = "50.0 mm"

    # PMT outer layer repeater
    pmtOuterRepeaterType = "cubicArray"
    pmtOuterRepeaterX = "1"
    pmtOuterRepeaterY = "2"
    pmtOuterRepeaterZ = "2"
    pmtOuterRepeaterVector = "0.0 18.04 18.94 mm"

    # PMT inner layer
    pmtInnerMaterial = "Vacuum"
    pmtInnterTranslation = "0.0 0.0 0.0 mm"
    pmtInnerRadiusMax = "7.79 mm"
    pmtInnerRadiusMin = "0.00 mm"
    pmtInnerHeight = "48.0 mm"

    # PMT top cap
    pmtCapMaterial = "Borosilicate-Glass"
    pmtTCTranslation = "0.0 0.0 24.5 mm"
    pmtBCTranslation = "0.0 0.0 -24.5 mm"
    pmtCapRadiusMax = "7.79 mm"
    pmtCapRadiusMin = "0.00 mm"
    pmtCapHeight = "1.0 mm"

    # Phantom name
    phantomName1 = "phantom1"
    activatePhantom2 = False
    phantomName2 = "phantom2"

    # cylinder phantom 1
    cylinderPhantom1 = False
    cylinderPhantomRadiusMin1 = "0.0 cm"
    cylinderPhantomRadiusMax1 = "5.0 mm"
    cylinderPhantomHeight1 = "11.5 mm"
    cylinderPhantomMaterial1 = "PMMA"
    cylinderPhantomTranslation1 = "-45.83 55.10 -4.15 mm" # 4.75

    # NRW-100 phantom 1 
    nrwPhantom1 = False
    nrwPhantomRadiusMax1 = "{} um".format(tracerSize)
    nrwPhantomRadiusMin1 = "0.0 um"
    nrwPhantomMaterial1 = "Nrw-mat"
    nrwPhantomTranslation1 = "-400.0 -400.0 -85.0 mm"

    # More complex NRW-100 tracer/phantom
    nrwPhantomShell1 = True
    nrwPhantomShell1FullRmin = "0.0 um"
    nrwPhantomShell1FullRmax = "{} um".format(tracerSize + 85)
    nrwPhantomShell1Material = "CyanoAcrylate"
    nrwPhantomShell1Translation = "20.0 20.0 20.0 mm" # Tracer location.  NB VERY IMPORTANT !!!.  This doesn't matter when you input a placements file
    # Don't touch below 
    #With Above
    nrwPhantomSilicaRmin = "{} um".format(tracerSize)
    nrwPhantomSilicaRmax = "{} um".format(tracerSize + 25)
    nrwPhantomSilicaMaterial = "Quartz"
    # With Above
    nrwPhantomTracerRmin = " 0.0 um"
    nrwPhantomTracerRmax = "{} um".format(tracerSize)
    nrwPhantomTracerMaterial = "Nrw-mat"
    #With Above
    nrwPhantomTracerForbidRmin = "0.0 um"
    nrwPhantomTracerForbidRmax = "{} um".format(tracerSize - 10)
    nrwPhantomTracerForbidMaterial = "Nrw-mat"

    # cylinder phantom 2
    cylinderPhantom2 = False
    cylinderPhantomRadiusMin2 = "0.0 cm"
    cylinderPhantomRadiusMax2 = "5.0 mm"
    cylinderPhantomHeight2 = "11.5 mm"
    cylinderPhantomMaterial2 = "PMMA"
    cylinderPhantomTranslation2 = "0.00 70.00 0.00 mm"

    # NRW-100 phantom 2
    nrwPhantom2 = False
    nrwPhantomRadiusMax2 = "50.0 um"
    nrwPhantomRadiusMin2 = "0.0 um"
    nrwPhantomMaterial2 = "Nrw-mat"
    nrwPhantomTranslation2 = "100.00 100.00 40.00 mm"

    # STL phantom
    stlPhantom = False
    stlPhantomTranslation = "0.00 0.00 0.00 mm"
    stlPhantomPath = "./data/smallTubeCone.stl"
    stlPhantomMaterial = "PMMA"
    stlPhantomColour = "red"

    # Voxel phantom
    voxelPhantom = False
    voxelPhantomType = "ImageNestedParametrisedVolume"
    voxelPhantomImageLocation = "data/bubbleTest.hdr"
    voxelPhantomMaterial = "data/AttenuationRange.dat"
    voxelPhantomTranslation = "0.0 0.0 0.0 mm"
    voxelPhantomRotationAxis = "1 0 0"
    voxelPhantomRotationAngle = "0 deg"

    # cylinderMaterial = "{}".format(Media[iv])
    # cylinderRadMax = "300.0 mm"
    # cylinderRadMin = "0.0 mm"
    # CylinderHeight = "230.0 mm"

    # Moving Phantom Controls
    movePhantom1 = True
    movePhantom2 = False

    #   MOTION 1
    # Generic Motion
    genericMotion1 = True
    genericMotionFileName1 = "/home/rayhaan/REPO_HR++/GATE_HR/Circles/Placements/cirlce_vT_10.0mps_t_2period_deltaT_1e-05s_radius_{}mm.placements".format(radii[h])
    # linear Motion
    linearMotion1 = False
    linearMotionSpeed1 = "400.0 400.0 85.0 mm/s"

    # orbiting motion
    orbitingMotion1 = False
    orbitingMotionSpeed1 = "382.278222 deg/s"
    orbitingMotionSetPoint11 = "-1.30 -0.54 -1.0 mm"
    orbitingMotionSetPoint21 = "-1.30 -0.54 1.0 mm"

    # oscilating motion
    oscilatingMotion1 = False
    oscilatingMotionAmplitude1 = "0.0 0.0 0.54 mm"
    oscilatingMotionFrequency1 = "1.06 Hz"
    oscilatingMotionSetPhase1 = "-1.28 rad" # 2.18

    #   MOTION 2
    # Generic Motion
    genericMotion2 = False
    genericMotionFileName2 = "PresentationTurbulentPath.placements"

    # linear Motion
    linearMotion2 = False
    linearMotionSpeed2 = "-400.0 -400.0 -85.0 mm/s"

    # orbiting motion
    orbitingMotion2 = False
    orbitingMotionSpeed2 = "1800.0  deg/s"
    orbitingMotionSetPoint12 = "0.0 0.0 -1.0 mm"
    orbitingMotionSetPoint22 = "0.0 0.0 1.0 mm"

    # oscilating motion
    oscilatingMotion2 = False
    oscilatingMotionAmplitude2 = "0.0 0.44211748 0.0 mm"
    oscilatingMotionFrequency2 = "0.42317 Hz"
    oscilatingMotionSetPhase2 = "3.14159265359 rad"

    #   PHYSICS AND PRODUCTION CUTS
    # crystal cuts
    crystalGammaCut = "crystal 1.0 cm"
    crystalElectronCut = "crystal 1.0 cm"
    crystalPositronCut = "crystal 1.0 cm"

    #   PHANTOM CUTS
    # phantom cuts 1
    phantomGammaCut1 = "phantom1 0.1 mm"
    phantomElectronCut1 = "phantom1 0.1 mm"
    phantomPositronCut1 = "phantom1 0.1 mm"
    phantomMaxStepSize1 = "phantom1 0.01 mm"

    phantomGammaCut2 = "phantom2 0.1 mm"
    phantomElectronCut2 = "phantom2 0.1 mm"
    phantomPositronCut2 = "phantom2 0.1 mm"
    phantomMaxStepSize2 = "phantom2 0.01 mm"

    #   DIGITISER
    # readout
    readoutDepth = "1"
    readoutPolicy = "TakeEnergyCentroid"

    # energy blurring
    energyBlurringMinRes = "0.20"
    energyBlurringMaxRes = "0.30"
    energyBlurringQE = "0.88"
    energyBlurringReferenceEnergy = "511.0 keV"

    # energy cuts
    energyCutThreshold = "350.0 keV"
    energyCutUphold = "650.0 keV"

    # timing resolution
    timingResolution = "2.0 ns"

    # noise
    NOISE = False
    NoiseEnergyMean = "511 keV"
    NoiseEnergySigma = "0.35 keV"
    NoiseLambda = "1.85 us" # 0.85

    # deadtime
    deadtimeValue = "700.0 ns"
    deadtimeMode = "paralysable"
    deadtimeVolume = "block"

    # coincidence sorter
    coincidenceWindow = "12. ns"
    coincidenceWindowOffset = "0.0 ns"
    coincidenceWindowMineSectorDifference = "19"
    coincidenceWindowPolicy = "keepIfAnyIsGood"

    # delayed coincidence sorter
    delayedWindow = "64000.0 ns"
    delayedWindowOffset = "12.0 ns"

    #   SOURCES
    # Source 1
    ga68source1 = False
    ga68sourceSphere1 = True # For the complex NRW-100 tracer
    na22source1 = False
    f18source1 = False
    sourceActivity1 = "32967000 Bq" # Activity of tracer.  Correct for electron capture!!!  
    sourceAttach1 = "nrw100Tracer"
    sourceForbid1 = "nrw100TracerForbid"
    sourceActivity1Radius = "{} um".format(tracerSize)

    # Source 2
    ga68source2 = False
    na22source2 = False
    f18source2 = False
    sourceActivity2 = "16483500 Bq"
    sourceAttach2 = "phantom2"

    # Voxel Source, don't bother with this if you can help it.  A bit complex to use.  
    voxelSource = False
    voxelSourcePixel2Activity = "data/ActivityRange.dat"
    voxelSourceImage = "data/bubbleTest.hdr"
    voxelSourcePosition = "-128.0 -128.0 -100.0 mm"
    voxelSourceType = "backtoback"
    voxelSourceParticle = "gamma"
    voxelSourceParticleEnergyType = "Mono"
    voxelSourceParticleEnergy = "511.0 keV"
    voxelSourceAngularDist = "iso"

    #   TIMING CONTROLS
    startTime = 0.0 # Simulation start time in s
    endTime = totalTime # Sim end time in s
    timeStep = "0.000001 s" # Make this less than your custom path time step or the desired temporal resolution if using the predefined motions
    coreTime = endTime/cores
    coreCounter = 0


    ############################
    ###### CREATE MACRO/S ######
    ############################

    for i in range(cores):
        coreCounterStart = coreCounter
        coreCounterEnd = coreCounter + coreTime
        #Create and open file
        file = open("{}{}_{}.mac".format(MACRO_FOLDER,NAME, i),"w")
        

        
        file.write(
            "#  WORLD\n"
            "/gate/geometry/setMaterialDatabase                      data/GateMaterials.db                # Database of materials that can be used\n"
            "/gate/world/geometry/setXLength                         {}                                   # The x component of the world\n"
            "/gate/world/geometry/setYLength                         {}                                   # The y component of the world\n"
            "/gate/world/geometry/setZLength                         {}                                   # The z component of the world\n"
            "/gate/world/setMaterial                                 {}                                   # The material the world\n"
            "\n".format(worldXLength, worldYLength, worldZLength, worldMaterial)
            )

        file.write(
            "#   CPET\n"
            "/gate/world/daughters/name                              CPET                                 # The name of the new volume\n" 
            "/gate/world/daughters/insert                            cylinder                             # Inserting a cylinder into the world volume\n"
            "/gate/CPET/setMaterial                                  {}                                   # The material the cylinder is made of\n" 
            "/gate/CPET/geometry/setRmax                             {}                                   # The outer radius of the cylinder\n"
            "/gate/CPET/geometry/setRmin                             {}                                   # The inner radius of the cylinder\n"
            "/gate/CPET/geometry/setHeight                           {}                                   # The height/length of the cylinder\n"
            "/gate/CPET/vis/forceWireframe                                                                # Forces the volume to be a wireframe.  It is transparent\n"
            "\n".format(CPETcylinderMaterial, CPETcylinderMaxRadius, CPETcylinderMinRadius, CPETclyinderHeight)
        )

        file.write(
            "#   Sector\n"
            "/gate/CPET/daughters/name                              sector                                # The name of the new volume\n"
            "/gate/CPET/daughters/insert                            cylinder                              # Inserting a cylinder into the world volume\n"
            "/gate/sector/setMaterial                             {}                                    # The material the cylinder is made of\n"
            "/gate/sector/placement/setTranslation                {}                          # The position of the box\n"
            "/gate/sector/geometry/setRmax                        {}                                # The outer radius of the cylinder\n" 
            "/gate/sector/geometry/setRmin                        {}                                # The inner radius of the cylinder\n"
            "/gate/sector/geometry/setHeight                      {}                                 # The height/length of the cylinder\n"
            "/gate/sector/geometry/setPhiStart                    {}                                # The start angle of the cylinder\n"
            "/gate/sector/geometry/setDeltaPhi                    {}                                   # The phi or x-y angle that the cylinder will span starting from phi start and ending at phi star + delta phi\n"
            "/gate/sector/vis/forceWireframe                                                              # Forces the volume to be a wireframe.  It is transparent\n"
            "/gate/sector/vis/setVisible                          1                                       # Visibility option\n"
            "/gate/sector/vis/setColor                            red                                     # Setting the colour of the geometry\n"
            "\n".format(sectorMaterial, sectorTranslation, sectorMaxRadius, sectorMinRadius, sectorHeight, sectorPhiStart, sectorDeltaPhi)
        )

        # file.write(
        #     "#  SECTOR REPEAT\n"
        #     "/gate/sector/repeaters/insert                         {}                                 # Linear repeater.  Repeats geometry in one direction with centre-to-centre spacing\n"
        #     "/gate/sector/ring/setRepeatNumber                     {}                                    # The number of times the geometry is going to be repeated\n"
        #     "/gate/sector/ring/excludeBlock1                       {}                                     # The fist blocks segment to be excluded, a max of 4 can be removed\n" 
        #     "/gate/sector/ring/excludeBlock2                       {}                                    # The second blocks segment to be excluded, a max of 4 can be removed\n"
        #     "\n".format(sectorRepeatType, sectorRepeatNo, sectorRepeatExclude1, sectorRepeatExclude2)
        # )

        file.write(
            "#  SECTOR REPEAT\n"
            "/gate/sector/repeaters/insert                         {}                                 # Linear repeater.  Repeats geometry in one direction with centre-to-centre spacing\n"
            "/gate/sector/ring/setRepeatNumber                     {}                                    # The number of times the geometry is going to be repeated\n"
            # "/gate/sector/ring/excludeBlock1                       {}                                     # The fist blocks segment to be excluded, a max of 4 can be removed\n" 
            # "/gate/sector/ring/excludeBlock2                       {}                                    # The second blocks segment to be excluded, a max of 4 can be removed\n"
            "\n".format(sectorRepeatType, sectorRepeatNo)
        )

        file.write(
            "#   CASSETTE\n"
            "/gate/sector/daughters/name                            cassette                             # The name of the new volume\n"
            "/gate/sector/daughters/insert                          cylinder                             # Inserting a cylinder into the world volume\n"
            "/gate/cassette/setMaterial                             {}                                  # The material the cylinder is made of\n"
            "/gate/cassette/placement/setTranslation                {}                      # The position of the box\n"
            "/gate/cassette/geometry/setRmax                        {}                             # The outer radius of the cylinder\n"
            "/gate/cassette/geometry/setRmin                        {}                             # The inner radius of the cylinder\n"
            "/gate/cassette/geometry/setHeight                      {}                              # The height/length of the cylinder\n"
            "/gate/cassette/geometry/setPhiStart                    {}                             # ++ Same as Above ++\n"
            "/gate/cassette/geometry/setDeltaPhi                    {}                                # ++ Same as Above ++\n"
            "/gate/cassette/vis/forceWireframe                                                           # Forces the volume to be a wireframe.  It is transparent\n"
            "/gate/cassette/vis/setColor                            cyan                                 # ++ Same as Above ++\n"
            "/gate/cassette/vis/setVisible                          0                                    # ++ Same as Above ++\n"
            "\n".format(cassetteMaterial, cassetteTranslation, cassetteMaxRadius, cassetteMinRadius, cassetteHeight, cassettePhiStart, cassetteDeltaPhi)
        )

        file.write(
            "#   BLOCK\n"
            "/gate/cassette/daughters/name                           block                               # ++ Same as Above ++\n"
            "/gate/cassette/daughters/insert                         box                                 # Defining a box that exists within ecat\n"
            "/gate/block/placement/setTranslation                    {}                                  # The position of the box 449.85 -21.37449981\n"
            "/gate/block/placement/setRotationAxis                   {}                               # ++ Same as Above ++\n"
            "/gate/block/placement/setRotationAngle                  {}                            # ++ Same as Above ++\n"
            "/gate/block/geometry/setXLength                         {}                             # The x component of the box\n"
            "/gate/block/geometry/setYLength                         {}                            # The y component of the box\n"
            "/gate/block/geometry/setZLength                         {}                            # The z component of the box\n"
            "/gate/block/setMaterial                                 {}                                 # The material which the box is made of\n"
            "/gate/block/vis/forceWireframe                                                              # Forces the volume to be a wireframe.  It is transparent\n"
            "/gate/block/vis/setVisible                              1\n"
            "\n".format(blockTranslation, blockRotationAxis, blockRotationAngle, blockXLength, blockYLength, blockZLength, blockMaterial)
        )

        file.write(
            "#   BLOCK REPEAT\n"
            "/gate/block/repeaters/insert                            {}                              # Linear repeater.  Repeats geometry in one direction with centre-to-centre spacing\n"
            "/gate/block/linear/setRepeatNumber                      {}                                  # The number of times the geometry is going to be repeated\n"
            "/gate/block/linear/setRepeatVector                      {}                     # The repeat vector specifies the vector centre-to-centre spacing of the repeated volumes\n"
            "\n".format(blockRepeatType, blockRepeatNo, blockRepeatVector)
        )

        file.write(
            "#   CRYSTAL\n"
            "/gate/block/daughters/name                              crystal                             # ++ Same as Above ++\n"
            "/gate/block/daughters/insert                            box                                 # ++ Same as Above ++\n"
            "/gate/crystal/geometry/setXLength                       {}                             # ++ Same as Above ++\n"
            "/gate/crystal/geometry/setYLength                       {}                             # ++ Same as Above ++\n"
            "/gate/crystal/geometry/setZLength                       {}                             # ++ Same as Above ++\n"
            "/gate/crystal/placement/setTranslation                  {}                    # ++ Same as Above ++\n"
            "/gate/crystal/setMaterial                               {}                                 # ++ Same as Above ++\n"
            "/gate/crystal/vis/setColor                              yellow                              # ++ Same as Above ++\n"
            "\n".format(crystalXLength, crystalYLength, crystalZLength, crystalTranslation, crystalMaterial)
        )

        file.write(
            "#   REPEAT CRYSTAL\n"
            "/gate/crystal/repeaters/insert                          {}                          # A cubic array repeater, Repeats volumes in the x,y and z direction\n"              
            "/gate/crystal/cubicArray/setRepeatNumberX               {}                                   # The number of x repitiions\n" 
            "/gate/crystal/cubicArray/setRepeatNumberY               {}                                   # The number of y repitiions\n"
            "/gate/crystal/cubicArray/setRepeatNumberZ               {}                                   # The number of z repitiions\n"
            "/gate/crystal/cubicArray/setRepeatVector                {}                     # The repeat vector specifies the centre-to-centre spacing of the repeated volumes,  0. 4.51 4.85\n"
            "\n".format(crystalRepeatType, crystalRepeatX, crystalRepeatY, crystalRepeatZ, crystalRepeatVector)
        )

        file.write(
            "#   PMT OUTER LAYER\n"
            "/gate/block/daughters/name                             pmt1                                 # ++ Same as Above ++\n" 
            "/gate/block/daughters/insert                           cylinder                             # ++ Same as Above ++\n"
            "/gate/pmt1/placement/setTranslation                    {}                      # ++ Same as Above ++\n" 
            "/gate/pmt1/placement/alignToX                                                               # Aligns the axial component of the cylinder with the x axis\n"
            "/gate/pmt1/geometry/setRmax                            {}                              # ++ Same as Above ++\n"
            "/gate/pmt1/geometry/setRmin                            {}                              # ++ Same as Above ++\n"
            "/gate/pmt1/geometry/setHeight                          {}                              # ++ Same as Above ++\n"
            "/gate/pmt1/setMaterial                                 {}                   # ++ Same as Above ++\n"
            "/gate/pmt1/vis/forceSolid                                                                   # Forces the volume to be a solid and not wireframe\n"
            "/gate/pmt1/vis/setVisible                              1                                    # ++ Same as Above ++\n"
            "/gate/pmt1/vis/setColor                                red                                  # ++ Same as Above ++\n"
            "\n".format(pmtOuterTranslation, pmtOuterRadiusMax, pmtOuterRadiusMin, pmtOuterHeight, pmtOuterMaterial)
        )

        file.write(
            "#	PMT OUTER LAYER REPEATER\n"
            "/gate/pmt1/repeaters/insert                            {}                           # ++ Same as Above ++\n"
            "/gate/pmt1/cubicArray/setRepeatNumberX                 {}                                    # ++ Same as Above ++\n"
            "/gate/pmt1/cubicArray/setRepeatNumberY                 {}                                    # ++ Same as Above ++\n"
            "/gate/pmt1/cubicArray/setRepeatNumberZ                 {}                                    # ++ Same as Above ++\n"  
            "/gate/pmt1/cubicArray/setRepeatVector                  {}                   # ++ Same as Above ++\n"
            "\n".format(pmtOuterRepeaterType, pmtOuterRepeaterX, pmtOuterRepeaterY, pmtOuterRepeaterZ, pmtOuterRepeaterVector)
        )

        file.write(
            "#   PMT INNER LAYER\n"
            "/gate/pmt1/daughters/name                              pmt2                                 # ++ Same as Above ++\n"
            "/gate/pmt1/daughters/insert                            cylinder                             # ++ Same as Above ++\n"
            "/gate/pmt2/placement/setTranslation                    {}                       # ++ Same as Above ++\n"
            "#/gate/pmt2/placement/alignToX                                                              # ++ Same as Above ++\n"
            "/gate/pmt2/geometry/setRmax                            {}                              # ++ Same as Above ++\n"
            "/gate/pmt2/geometry/setRmin                            {}                              # ++ Same as Above ++\n"
            "/gate/pmt2/geometry/setHeight                          {}                              # ++ Same as Above ++\n"
            "/gate/pmt2/setMaterial                                 {}                               # ++ Same as Above ++\n"
            "/gate/pmt2/vis/forceSolid                                                                   # ++ Same as Above ++\n"
            "/gate/pmt2/vis/setVisible                              0                                    # ++ Same as Above ++\n"
            "/gate/pmt2/vis/setColor                                green                                # ++ Same as Above ++\n"
            "\n".format(pmtInnterTranslation, pmtInnerRadiusMax, pmtInnerRadiusMin, pmtInnerHeight, pmtInnerMaterial)
        )

        file.write(
            "#   PMT TOP CAP\n"
            "/gate/pmt1/daughters/name                              pmt1TopCap                           # ++ Same as Above ++\n"
            "/gate/pmt1/daughters/insert                            cylinder                             # ++ Same as Above ++\n"
            "/gate/pmt1TopCap/placement/setTranslation              {}                      # ++ Same as Above ++\n"
            "#/gate/pmt2/placement/alignToX                                                              # ++ Same as Above ++\n"
            "/gate/pmt1TopCap/geometry/setRmax                      {}                             # ++ Same as Above ++\n"
            "/gate/pmt1TopCap/geometry/setRmin                      {}                              # ++ Same as Above ++\n"
            "/gate/pmt1TopCap/geometry/setHeight                    {}                              # ++ Same as Above ++\n"
            "/gate/pmt1TopCap/setMaterial                           {}                    # ++ Same as Above ++\n"
            "/gate/pmt1TopCap/vis/forceSolid                                                             # ++ Same as Above ++\n"
            "/gate/pmt1TopCap/vis/setVisible                        1                                    # ++ Same as Above ++\n"
            "/gate/pmt1TopCap/vis/setColor                          red                                  # ++ Same as Above ++\n"
            "\n".format(pmtTCTranslation, pmtCapRadiusMax, pmtCapRadiusMin, pmtCapHeight, pmtCapMaterial)
        )

        file.write(
            "#   PMT BOTTOM CAP\n"
            "/gate/pmt1/daughters/name                              pmt1BottomCap                        # ++ Same as Above ++\n"
            "/gate/pmt1/daughters/insert                            cylinder                             # ++ Same as Above ++\n"
            "/gate/pmt1BottomCap/placement/setTranslation           {}                     # ++ Same as Above ++\n"
            "#/gate/pmt2/placement/alignToX                                                              # ++ Same as Above ++\n"
            "/gate/pmt1BottomCap/geometry/setRmax                   {}                              # ++ Same as Above ++\n"
            "/gate/pmt1BottomCap/geometry/setRmin                   {}                              # ++ Same as Above ++\n"
            "/gate/pmt1BottomCap/geometry/setHeight                 {}                               # ++ Same as Above ++\n"
            "/gate/pmt1BottomCap/setMaterial                        {}                         # ++ Same as Above ++\n"
            "/gate/pmt1BottomCap/vis/forceSolid                                                          # ++ Same as Above ++\n"
            "/gate/pmt1BottomCap/vis/setVisible                     1                                    # ++ Same as Above ++\n"
            "/gate/pmt1BottomCap/vis/setColor                       red                                  # ++ Same as Above ++\n"
            "\n".format(pmtBCTranslation, pmtCapRadiusMax, pmtCapRadiusMin, pmtCapHeight, pmtCapMaterial)
        )

        file.write(
            "#   ATTACH SYSTEM\n"
            "/gate/systems/CPET/sector/attach                       sector                               # GATE has a template for detector systems.  CPET is one.  The components are ordered in a certain way.  For exmaple crystals and blocks need to be ordered in a certain manner\n"
            "/gate/systems/CPET/cassette/attach                     cassette                             # Adding the cassette component of the CPET system\n"
            "/gate/systems/CPET/module/attach                       block                                # Adding the block component of the CPET system\n"
            "/gate/systems/CPET/crystal/attach                      crystal                              # Adding the crystal component of the CPET system.\n"
            "\n"
        )

        file.write(
            "#   ATTACH CRYSTAL SD\n"
            "/gate/crystal/attachCrystalSD                                                               # Making the Crystal Sensitive to interactions\n"
            "\n"
        )

        if cylinderPhantom1 == True:
            file.write(
                "#   PHANTOM 1\n"
                "/gate/world/daughters/name         		               phantom1                              # ++ Same as Above ++\n" 
                "/gate/world/daughters/insert       		               cylinder                             # ++ Same as Above ++\n" 
                "/gate/phantom1/geometry/setRmin     		               {}                               # ++ Same as Above ++\n" 
                "/gate/phantom1/geometry/setRmax     		               {}                               # ++ Same as Above ++\n" 
                "/gate/phantom1/geometry/setHeight   	                   {}                              # ++ Same as Above ++\n" 
                "/gate/phantom1/setMaterial          		               {}                                 # ++ Same as Above ++\n" 
                "/gate/phantom1/vis/setColor         		               red                                  # ++ Same as Above ++\n" 
                "/gate/phantom1/placement/setTranslation                    {}                   # Specify the placement of the source within the world volume\n"
                "\n".format(cylinderPhantomRadiusMin1, cylinderPhantomRadiusMax1, cylinderPhantomHeight1, cylinderPhantomMaterial1, cylinderPhantomTranslation1)
            )

        if nrwPhantom1 == True:
            file.write(
                "#   NRW-100 PHANTOM 1\n"
                "/gate/world/daughters/name                                 phantom1\n"
                "/gate/world/daughters/insert                               sphere\n"
                "/gate/phantom1/geometry/setRmin                             {}\n"
                "/gate/phantom1/geometry/setRmax                             {}\n"
                "/gate/phantom1/setMaterial                                  {}\n"
                "/gate/phantom1/vis/setColor                                 red\n"
                "/gate/phantom1/placement/setTranslation                     {}\n"
                "\n".format(nrwPhantomRadiusMin1, nrwPhantomRadiusMax1, nrwPhantomMaterial1, nrwPhantomTranslation1)   
            )

        if nrwPhantomShell1 == True:
            file.write(
                "#  ACRYLATE REGION PHANTOM 1\n"
                "/gate/world/daughters/name                                 phantom1\n"
                "/gate/world/daughters/insert                                sphere\n"
                "/gate/phantom1/geometry/setRmin                             {}\n"
                "/gate/phantom1/geometry/setRmax                             {}\n"
                "/gate/phantom1/setMaterial                                  {}\n"
                "/gate/phantom1/vis/setColor                                 red\n"
                "/gate/phantom1/placement/setTranslation                     {}\n"
                "\n".format(nrwPhantomShell1FullRmin, nrwPhantomShell1FullRmax, nrwPhantomShell1Material, nrwPhantomShell1Translation)
            )

            file.write(
                "#  SILICA PHANTOM 1\n"
                "/gate/phantom1/daughters/name                            silicaRegion\n"
                "/gate/phantom1/daughters/insert                                sphere\n"
                "/gate/silicaRegion/geometry/setRmin                             {}\n"
                "/gate/silicaRegion/geometry/setRmax                             {}\n"
                "/gate/silicaRegion/setMaterial                                  {}\n"
                "/gate/silicaRegion/vis/setColor                                 red\n"
                "/gate/silicaRegion/placement/setTranslation                     0.0 0.0 0.0 mm\n"
                "\n".format(nrwPhantomSilicaRmin, nrwPhantomSilicaRmax, nrwPhantomSilicaMaterial)
            )

            file.write(
                "#   NRW-100 PHANTOM 1\n"
                "/gate/phantom1/daughters/name                                 nrw100Tracer\n"
                "/gate/phantom1/daughters/insert                               sphere\n"
                "/gate/nrw100Tracer/geometry/setRmin                             {}\n"
                "/gate/nrw100Tracer/geometry/setRmax                             {}\n"
                "/gate/nrw100Tracer/setMaterial                                  {}\n"
                "/gate/nrw100Tracer/vis/setColor                                 red\n"
                "/gate/nrw100Tracer/placement/setTranslation                     0.0 0.0 0.0 mm\n"
                "\n".format(nrwPhantomTracerRmin, nrwPhantomTracerRmax, nrwPhantomTracerMaterial)
            )

            file.write(
                "#   NRW-100 PHANTOM - FORBID 1\n"
                "/gate/nrw100Tracer/daughters/name                                 nrw100TracerForbid\n"
                "/gate/nrw100Tracer/daughters/insert                               sphere\n"
                "/gate/nrw100TracerForbid/geometry/setRmin                             {}\n"
                "/gate/nrw100TracerForbid/geometry/setRmax                             {}\n"
                "/gate/nrw100TracerForbid/setMaterial                                  {}\n"
                "/gate/nrw100TracerForbid/vis/setColor                                 red\n"
                "/gate/nrw100TracerForbid/placement/setTranslation                     0.0 0.0 0.0 mm\n"
                "\n".format(nrwPhantomTracerForbidRmin, nrwPhantomTracerForbidRmax, nrwPhantomTracerForbidMaterial)
            )
        
        if stlPhantom == True:
            file.write(
                "#   STL - PHANTOM\n"
                "/gate/world/daughters/name                                  phantom1\n"
                "/gate/world/daughters/insert                                tessellated\n"
                "/gate/phantom1/placement/setTranslation                      {}\n"
                "/gate/phantom1/geometry/setPathToSTLFile                     {}\n"
                "/gate/phantom1/setMaterial                                   {}\n"
                "/gate/phantom1/vis/setColor                                  {}\n"
                "/gate/phantom1/vis/forceWireframe\n"
                "\n".format(stlPhantomTranslation, stlPhantomPath, stlPhantomMaterial, stlPhantomColour)
            )
        
        if activatePhantom2 == True:
            if cylinderPhantom2 == True:
                file.write(
                    "#   PHANTOM 2\n"
                    "/gate/world/daughters/name         		               phantom2                              # ++ Same as Above ++\n" 
                    "/gate/world/daughters/insert       		               cylinder                             # ++ Same as Above ++\n" 
                    "/gate/phantom2/geometry/setRmin     		               {}                               # ++ Same as Above ++\n" 
                    "/gate/phantom2/geometry/setRmax     		               {}                               # ++ Same as Above ++\n" 
                    "/gate/phantom2/geometry/setHeight   	                   {}                              # ++ Same as Above ++\n" 
                    "/gate/phantom2/setMaterial          		               {}                                 # ++ Same as Above ++\n" 
                    "/gate/phantom2/vis/setColor         		               red                                  # ++ Same as Above ++\n" 
                    "/gate/phantom2/placement/setTranslation                    {}                   # Specify the placement of the source within the world volume\n"
                    "\n".format(cylinderPhantomRadiusMin2, cylinderPhantomRadiusMax2, cylinderPhantomHeight2, cylinderPhantomMaterial2, cylinderPhantomTranslation2)
                )

            if nrwPhantom2 == True:
                file.write(
                    "#   NRW-100 PHANTOM 2\n"
                    "/gate/world/daughters/name                                 phantom2\n"
                    "/gate/world/daughters/insert                               sphere\n"
                    "/gate/phantom2/geometry/setRmin                             {}\n"
                    "/gate/phantom2/geometry/setRmax                             {}\n"
                    "/gate/phantom2/setMaterial                                  {}\n"
                    "/gate/phantom2/vis/setColor                                 red\n"
                    "/gate/phantom2/placement/setTranslation                     {}\n"
                    "\n".format(nrwPhantomRadiusMin2, nrwPhantomRadiusMax2, nrwPhantomMaterial2, nrwPhantomTranslation2)   
                )



        
        if voxelPhantom == True:
            file.write(
                "/gate/world/daughters/name                                  phantom1\n"
                "/gate/world/daughters/insert                                {}\n"
                "/gate/phantom1/geometry/setImage                             {}\n"
                "/gate/phantom1/geometry/setRangeToMaterialFile               {}\n"
                "/gate/phantom1/placement/setTranslation                      {}\n"
                "/gate/phantom1/placement/setRotationAxis                     {}\n"
                "/gate/phantom1/placement/setRotationAngle                    {}\n"
                "\n".format(voxelPhantomType, voxelPhantomImageLocation, voxelPhantomMaterial, voxelPhantomTranslation, voxelPhantomRotationAxis, voxelPhantomRotationAngle)
            )
        
        if movePhantom1 == True:
            if genericMotion1 == True:
                file.write(
                    "#   GENERIC MOVE\n"
                    "/gate/phantom1/moves/insert                                 genericMove\n"
                    "/gate/phantom1/genericMove/setPlacementsFilename  {}\n"
                    "\n".format(genericMotionFileName1)
                )

            if linearMotion1 == True:
                file.write(
                    "# LINEAR MOTION\n"
                    "/gate/phantom1/moves/insert translation\n"
                    "/gate/phantom1/translation/setSpeed {}\n"
                    "\n".format(linearMotionSpeed1)
                )
            if orbitingMotion1 == True:
                file.write(
                    "#   ORBITING PHANTOM\n"
                    "/gate/phantom1/moves/insert                            orbiting                             # Specify the motion to simulate\n"
                    "/gate/phantom1/orbiting/setSpeed                       {}                        # Specify the speed of orbital motion -360 deg/s\n"
                    "/gate/phantom1/orbiting/setPoint1                      {}                    # Anchor point 1 for orbital axis -2.60415503 1.18758883 0.0\n"
                    "/gate/phantom1/orbiting/setPoint2                      {}                     # Anchor point 2 for orbital axis. -2.60415503 1.18758883 1.0\n"
                    "\n".format(orbitingMotionSpeed1, orbitingMotionSetPoint11, orbitingMotionSetPoint21)
                )
            
            if oscilatingMotion1 == True:
                file.write(
                    "#   OSCILATING PHANTOM\n"
                    "/gate/phantom1/moves/insert                            osc-trans                            # +++ Same as Above +++\n"
                    "/gate/phantom1/osc-trans/setAmplitude                  {}                # Set the amplitude of the oscillation\n"
                    "/gate/phantom1/osc-trans/setFrequency                  {}                           # Set the frequency of oscillation\n"
                    "/gate/phantom1/osc-trans/setPhase                      {}                       # Set the phase of the oscillation\n"
                    "/gate/geometry/rebuild\n"
                    "\n".format(oscilatingMotionAmplitude1, oscilatingMotionFrequency1, oscilatingMotionSetPhase1)
                )
        
        if movePhantom2 == True:
            if genericMotion2 == True:
                file.write(
                    "#   GENERIC MOVE\n"
                    "/gate/phantom2/moves/insert                                 genericMove\n"
                    "/gate/phantom2/genericMove/setPlacementsFilename  data/{}\n"
                    "\n".format(genericMotionFileName2)
                )

            if linearMotion2 == True:
                file.write(
                    "# LINEAR MOTION\n"
                    "/gate/phantom2/moves/insert translation\n"
                    "/gate/phantom2/translation/setSpeed {}\n"
                    "\n".format(linearMotionSpeed2)
                )
            if orbitingMotion2 == True:
                file.write(
                    "#   ORBITING PHANTOM\n"
                    "/gate/phantom2/moves/insert                            orbiting                             # Specify the motion to simulate\n"
                    "/gate/phantom2/orbiting/setSpeed                       {}                        # Specify the speed of orbital motion -360 deg/s\n"
                    "/gate/phantom2/orbiting/setPoint1                      {}                    # Anchor point 1 for orbital axis -2.60415503 1.18758883 0.0\n"
                    "/gate/phantom2/orbiting/setPoint2                      {}                     # Anchor point 2 for orbital axis. -2.60415503 1.18758883 1.0\n"
                    "\n".format(orbitingMotionSpeed2, orbitingMotionSetPoint12, orbitingMotionSetPoint22)
                )
            
            if oscilatingMotion2 == True:
                file.write(
                    "#   OSCILATING PHANTOM\n"
                    "/gate/phantom2/moves/insert                            osc-trans                            # +++ Same as Above +++\n"
                    "/gate/phantom2/osc-trans/setAmplitude                  {}                # Set the amplitude of the oscillation\n"
                    "/gate/phantom2/osc-trans/setFrequency                  {}                           # Set the frequency of oscillation\n"
                    "/gate/phantom2/osc-trans/setPhase                      {}                       # Set the phase of the oscillation\n"
                    "/gate/geometry/rebuild\n"
                    "\n".format(oscilatingMotionAmplitude2, oscilatingMotionFrequency2, oscilatingMotionSetPhase2)
                )
        file.write(
            "#  PHYSICS\n"
            "/gate/physics/addPhysicsList                               emstandard_opt4                         # em standard opt4 selected, best for E range we're intereseted in\n"
            "\n"
        )

        if activatePhantom2 == True:
            file.write(
                "/gate/physics/Gamma/SetCutInRegion                                      {}                          # Deposit gamma ray secondaries energy if the secondaries expected range is less than 0.1 mm in the phantom.\n"
                "/gate/physics/Electron/SetCutInRegion                                   {}                          # Deposit electron secondaries energy if the secondaries expected range is less than 0.1 mm in the phantom.\n"
                "/gate/physics/Positron/SetCutInRegion                                   {}                          # Deposit positron secondaries energy if the secondaries expected range is less than 0.1 mm in the phantom.\n"
                "\n"
                "/gate/physics/SetMaxStepSizeInRegion                                    {}                         # The maximum step size a particle can have, in the phantom region.\n"
                "\n".format(phantomGammaCut2, phantomElectronCut2, phantomPositronCut2, phantomMaxStepSize2)
            )

        file.write(
            "#   INITIALISE\n"
            "/gate/run/initialize                                                    # Initialize the simulaion\n"
            "\n"
        )

        file.write(
            "#	ADDER\n"
            "/gate/digitizer/Singles/insert                                          adder                       # The adder sums all the hits within crystals\n"
            "\n"
        )

        file.write(
            "#   READOUT\n"
            "/gate/digitizer/Singles/insert                                          readout                     # The readout component of the singles over which geometry to sum the energies of the crystal\n"
            "/gate/digitizer/Singles/readout/setDepth                                {}                           # The energy of an individual event is the sum of energies in a block.\n"
            "/gate/digitizer/Singles/readout/setPolicy                               {}          # The readout position is determined by a weighting the indices of each crystal pulse by energy deposited to get the enrgy centroid.  The position is the centre of the crystal whose crystal index was selected.\n"
            "\n".format(readoutDepth, readoutPolicy)
        )

        file.write(
            "#   ENERGY BLURRING\n"
            "/gate/digitizer/Singles/insert                                          crystalblurring             # Crystal blurring component.  The energy resolution blurring factor is randomly chosen between the lower and upper threshold for each crystal.\n"
            "/gate/digitizer/Singles/crystalblurring/setCrystalResolutionMin         {}                        # Lower threshold of blurring of Energy in crystal\n"
            "/gate/digitizer/Singles/crystalblurring/setCrystalResolutionMax         {}                        # Upper threshold of blurring of Energy in crystal\n"
            "/gate/digitizer/Singles/crystalblurring/setCrystalQE                    {}                       # Quantum efficiency of crystal, for CP data 0.545\n"
            "/gate/digitizer/Singles/crystalblurring/setCrystalEnergyOfReference     {}                   # Blurring energy reference\n"
            "\n".format(energyBlurringMinRes, energyBlurringMaxRes, energyBlurringQE, energyBlurringReferenceEnergy)
        )

        file.write(
            "#   ENERGY CUTS\n"
            "/gate/digitizer/Singles/insert                                          thresholder                 # The energy thresholder sets the energy range which the detector will accept singles\n"
            "/gate/digitizer/Singles/thresholder/setThreshold                        {}                    # The lower limit of the energy range\n"
            "/gate/digitizer/Singles/insert                                          upholder                    # The upholder allows you to define an upper limit for the energy that will be accepted\n"
            "/gate/digitizer/Singles/upholder/setUphold                              {}                    # The upper limit of the energy range\n"
            "\n".format(energyCutThreshold, energyCutUphold)
        )

        file.write(
            "#   TIMING RESOLUTION\n"
            "/gate/digitizer/Singles/insert                                          timeResolution              # Implementing the timiing resolution of the detector\n"
            "/gate/digitizer/Singles/timeResolution/setTimeResolution                {}                      # The value of the timing resolution, 2*tau\n"
            "\n".format(timingResolution)
        )
        if NOISE == True:
            file.write(
                "#   NOISE\n"
                "/gate/distributions/name                                                energy_distrib              # We added a distribution called energy_dist\n"
                "/gate/distributions/insert                                              Gaussian                    # The shape of the distribution is Gaussian\n"
                "/gate/distributions/energy_distrib/setMean                              {}                     # The mean energy of the of the Gaussian\n"
                "/gate/distributions/energy_distrib/setSigma                             {}                    # The standard deviation of the Gaussian\n"
                "/gate/distributions/name                                                dt_distrib                  # Adding another distribution called dt_distrib\n"
                "/gate/distributions/insert                                              Exponential                 # The shape of the distribtion is Exponential\n"
                "/gate/distributions/dt_distrib/setLambda                                {}                      # The lamdda value, A*e^(-lambda*t)\n"
                "/gate/digitizer/Singles/insert                                          noise                       # Inserting the noise component of the singles\n"
                "/gate/digitizer/Singles/noise/setDeltaTDistribution                     dt_distrib                  # Inserting the delta t distribution\n"
                "/gate/digitizer/Singles/noise/setEnergyDistribution                     energy_distrib              # Inserting the energy distribution\n"
                "\n".format(NoiseEnergyMean, NoiseEnergySigma, NoiseLambda)
            )

        file.write(
            "#   DEAD TIME\n"
            "/gate/digitizer/Singles/insert                                          deadtime                    # The deadtime implementation\n"
            "/gate/digitizer/Singles/deadtime/setDeadTime                            {}                     # The value of the deadtime\n"
            "/gate/digitizer/Singles/deadtime/setMode                                {}                 # The type of deadtime selected, there are two types, paralysable and non-paralysable\n" 
            "/gate/digitizer/Singles/deadtime/chooseDTVolume                         {}                       # The volume over which to apply the deadtime\n"
            "\n".format(deadtimeValue, deadtimeMode, deadtimeVolume)
        )

        file.write(
            "#   COINCIDENCE SORTER\n"
            "/gate/digitizer/Coincidences/setWindow                                  {}                      # The coincidence window implementation\n"
            "/gate/digitizer/Coincidences/setOffset                                  {}                       # The offset to the coincidence window\n"
            "/gate/digitizer/Coincidences/minSectorDifference                        {}                          # The minSectorDifference is similar to the span.  A minSectorDifference of 36, since we have 36 axial block segments, would be a prompt coincidence detector.\n" 
            "/gate/digitizer/Coincidences/MultiplesPolicy                            {}             # If at least one pair is good, keep the coincidence\n"
            "/gate/digitizer/Coincidences/describe                                                               # Print out for the coincidences\n"
            "\n".format(coincidenceWindow, coincidenceWindowOffset, coincidenceWindowMineSectorDifference, coincidenceWindowPolicy)
        )

        file.write(
            "#   DELAYED COINCIDENCE SORTER - Need the offset to work, choose 12 ns but might be incorrect\n"
            "/gate/digitizer/name                                                    delay                       # The name of the digitizer module\n"
            "/gate/digitizer/insert                                                  coincidenceSorter           # What type of module we are selecting\n" 
            "/gate/digitizer/delay/setWindow                                         {}                   # The value of the delayed coincidence window\n"
            "/gate/digitizer/delay/setOffset                                         {}                      # The delayed coincidence window offset.\n"  
            "/gate/digitizer/delay/describe                                                                      # Print out for the delays, should be 0\n"
            "\n".format(delayedWindow, delayedWindowOffset)
        )

        file.write(
            "#   SEED CONTROLLER\n"
            "/gate/random/setEngineName         JamesRandom                          # Method by which the seed will be generated\n"
            "#/gate/random/setEngineSeed        default                              # How you wish to set the engine seed, default\n"
            "#/gate/random/setEngineSeed        auto                                 # How you wish to set the engine seed, auto\n"
            "#/gate/random/setEngineSeed        123456789                            # Setting the engine seed manually\n"
            "/gate/random/setEngineSeed         default                              # ++ Same as Above ++\n"
            "#/gate/random/resetEngineFrom      fileName                             # Set the engine seed from a file\n"
            "/gate/random/verbose               1                                    # Set the random verbose level\n"
            "\n"
        )

        file.write(
            "#   STATS\n"
            "/gate/actor/addActor                                    SimulationStatisticActor stat           # Specify the actor to add and the name of the actor. The statistics actor outputs info about the simulation after every run.  Info like run #, event #, sim time, user time etc etc. See stat.txt for more info\n"
            "/gate/actor/stat/save                                   output/stat.txt                         # Specifies the output file of the statistics actor\n"
            "\n"
        )

        file.write(
            "#   ROOT OUTPUT\n"
            "/gate/output/root/enable                                                                        # Enable the root output.  Might add the *.npy output\n"
            "/gate/output/root/setFileName                           {}{}_{}               # The name and path of the root file CP1_3.4_1s_sim\n"
            "/gate/output/root/setRootSinglesAdderFlag               0                                       # Setting to include the SinglesAdder output as a TTREE\n"
            "/gate/output/root/setRootSinglesReadoutFlag             0                                       # Setting to include the SinglesReadout output as a TTREE\n"
            "/gate/output/root/setRootHitFlag                        0                                       # Setting to include the Hits output as a TTREE\n"
            "/gate/output/root/setRootSinglesFlag                    0                                       # Setting to include the Singles output as a TTREE\n"
            "/gate/output/root/setRootCoincidencesFlag               1                                       # Setting to include the Coincidences output as a TTREE\n"
            "\n".format(ROOT_FOLDER, "mp", i)
        )

        file.write(
            "#   ROOT SUMMARY\n"
            "/gate/output/summary/enable                                                                     # Enable summary simulation run\n"
            "/gate/output/summary/setFileName                        output/pet_summary.txt                  # Specify the ouput file of the summary\n"
            "/gate/output/summary/addCollection                      Singles                                 # The number of singles detected\n"
            "/gate/output/summary/addCollection                      Coincidences                            # The number of coincidences detected\n"
            "/gate/output/summary/addCollection                      delay                                   # The number of randoms detected\n"
            "\n"
        )

        file.write(
            "#   GATE VERBOSE OPTIONS\n"
            "/gate/verbose Physic    5                   # Physics verbose level\n"
            "/gate/verbose Cuts      0                   # production cuts verbose level\n"
            "/gate/verbose SD        0                   # Sensitive detector verbose level\n"
            "/gate/verbose Actions   0                   # Actions verbose level\n"
            "/gate/verbose Actor     0                   # Actor verbose level\n"
            "/gate/verbose Step      0                   # Step size verbose level\n"
            "/gate/verbose Error     0                   # Error verbose level\n"
            "/gate/verbose Warning   0                   # Warning verbose level\n"
            "/gate/verbose Output    0                   # Output verbose level\n"
            "/gate/verbose Beam      0                   # Beam verbose level\n"
            "/gate/verbose Volume    0                   # Volume verbose level\n"
            "/gate/verbose Image     0                   # Image verbose level\n"
            "/gate/verbose Geometry  0                   # Geometry verbose level\n"
            "/gate/verbose Core      0                   # Core verbose level\n"
            "\n"
            "#   GEANT4 VERBOSE OPTIONS\n"
            "/run/verbose            0                   # Run verbose level\n"
            "/event/verbose          0                   # Event verbose level\n"
            "/tracking/verbose       0                   # Tracking verbose level\n"
            "\n"
        )

        if ga68source1 == True:
            file.write(
                "#   GALLIUM - 68\n"
                "/gate/source/addSource                         ga68.1 gps                         # Source name\n"
                "/gate/source/ga68.1/setActivity                  {}                       # Activity of source  1978020 Bq\n"
                "/gate/source/ga68.1/gps/particle                 e+                               # Specifying the primary particle to be positrons\n"
                "/gate/source/ga68.1/gps/energytype               UserSpectrum                     # Specifying the energy type\n"
                "/gate/source/ga68.1/gps/setSpectrumFile          data/ga68_spectrum.txt           # Where to find user defined spectrum\n"
                "#/gate/source/ga68.1/gps/type                     Volume                          # Specifying that the source is a volume\n"
                "#/gate/source/ga68.1/gps/shape                    Sphere                          # Specifying the shape of the voxelised source\n"
                "#/gate/source/ga68.1/gps/radius                   301.0 um                         # Radius of source\n" 
                "#/gate/source/ga68.1/gps/radius0                  300.0 um                         # Radius of source\n" 
                "/gate/source/ga68.1/setForcedUnstableFlag        true                             # Forcing the source to be unstable and have the ability to decay\n"
                "/gate/source/ga68.1/setForcedHalfLife            4062.6 s                         # Sepecifying the half life of the source, 2.341e+7\n"
                "/gate/source/ga68.1/gps/centre                   0.0 0.0 0.0 cm                   # Specifying the position of the source\n"
                "/gate/source/ga68.1/gps/angtype                  iso                              # Specfying the angular distribution of the decay particles\n"
                "/gate/source/ga68.1/attachTo                     {}                          # Attaching the source to a volume.  This is for motion.  The location os the source is now relative to the volume it is attached to.\n"
                "\n".format(sourceActivity1, sourceAttach1)
            )

        if ga68sourceSphere1 == True:
            file.write(
                "#   GALLIUM - 68 - SPHERE\n"
                "/gate/source/addSource                         ga68_1 gps                         # Source name\n"
                "/gate/source/ga68_1/setActivity                  {}                       # Activity of source  1978020 Bq\n"
                "/gate/source/ga68_1/gps/particle                 e+                               # Specifying the primary particle to be positrons\n"
                "/gate/source/ga68_1/gps/energytype               UserSpectrum                     # Specifying the energy type\n"
                "/gate/source/ga68_1/gps/setSpectrumFile          data/ga68_spectrum.txt           # Where to find user defined spectrum\n"
                "/gate/source/ga68_1/gps/type                     Volume                          # Specifying that the source is a volume\n"
                "/gate/source/ga68_1/gps/shape                    Sphere                          # Specifying the shape of the voxelised source\n"
                "/gate/source/ga68_1/gps/radius                    {}                         # Radius of source\n" 
                "/gate/source/ga68_1/setForcedUnstableFlag        true                             # Forcing the source to be unstable and have the ability to decay\n"
                "/gate/source/ga68_1/setForcedHalfLife            4062.6 s                         # Sepecifying the half life of the source, 2.341e+7\n"
                "/gate/source/ga68_1/gps/centre                   0.0 0.0 0.0 cm                   # Specifying the position of the source\n"
                "/gate/source/ga68_1/gps/angtype                  iso                              # Specfying the angular distribution of the decay particles\n"
                "/gate/source/ga68_1/attachTo                     {}                          # Attaching the source to a volume.  This is for motion.  The location os the source is now relative to the volume it is attached to.\n"
                "/gate/source/ga68_1/gps/Forbid                        {}\n"
                "/gate/source/ga68_1/dump                          {}\n"
                "\n".format(sourceActivity1, sourceActivity1Radius ,sourceAttach1, sourceForbid1, 1)
            )


        if f18source1 == True:
            file.write(
                "#   FLUORINE - 18\n"
                "/gate/source/addSource                         f18.1 gps                         # ++ Same as Above ++\n"
                "/gate/source/f18.1/setActivity                   {}                      # ++ Same as Above ++\n" 
                "/gate/source/f18.1/gps/particle                  e+                              # ++ Same as Above ++\n" 
                "/gate/source/f18.1/gps/energytype                UserSpectrum                    # ++ Same as Above ++\n" 
                "/gate/source/f18.1/gps/setSpectrumFile           data/f18_spectrum.txt           # ++ Same as Above ++\n" 
                "/gate/source/f18.1/setForcedUnstableFlag         true                            # ++ Same as Above ++\n"  
                "/gate/source/f18.1/setForcedHalfLife             6586.26 s                       # ++ Same as Above ++\n"
                "/gate/source/f18.1/gps/centre                    0.0 0.0 0.0 cm                  # ++ Same as Above ++\n"
                "/gate/source/f18.1/gps/angtype                   iso                             # ++ Same as Above ++\n"
                "/gate/source/f18.1/attachTo                      {}                         # ++ Same as Above ++\n"
                "\n".format(sourceActivity1, sourceAttach1)
            )

        if na22source1 == True:
            file.write(
                "#   SODIUM - 22\n"
                "/gate/source/addSource                          na22.1 gps                        # ++ Same as Above ++\n"
                "/gate/source/na22.1/setActivity                   {}                      # ++ Same as Above ++\n"
                "/gate/source/na22.1/gps/particle                  e+                              # ++ Same as Above ++\n"
                "/gate/source/na22.1/gps/energytype                UserSpectrum                    # ++ Same as Above ++\n"
                "/gate/source/na22.1/gps/setSpectrumFile           data/na22_spectrum.txt          # ++ Same as Above ++\n"
                "/gate/source/na22.1/setForcedUnstableFlag         true                            # ++ Same as Above ++\n"
                "/gate/source/na22.1/setForcedHalfLife             8.199e+7 s                      # ++ Same as Above ++\n"
                "/gate/source/na22.1/gps/centre                    0.0 0.0 0.0 cm                  # ++ Same as Above ++\n"
                "/gate/source/na22.1/gps/angtype                   iso                             # ++ Same as Above ++\n"
                "/gate/source/na22.1/attachTo                      {}                         # ++ Same as Above ++\n"
                "\n".format(sourceActivity1, sourceAttach1)
            )

        if activatePhantom2 == True:
            if ga68source2 == True:
                file.write(
                    "#   GALLIUM - 68\n"
                    "/gate/source/addSource                         ga68.2 gps                         # Source name\n"
                    "/gate/source/ga68.2/setActivity                  {}                       # Activity of source  1978020 Bq\n"
                    "/gate/source/ga68.2/gps/particle                 e+                               # Specifying the primary particle to be positrons\n"
                    "/gate/source/ga68.2/gps/energytype               UserSpectrum                     # Specifying the energy type\n"
                    "/gate/source/ga68.2/gps/setSpectrumFile          data/ga68_spectrum.txt           # Where to find user defined spectrum\n"
                    "#/gate/source/ga68.2/gps/type                     Volume                          # Specifying that the source is a volume\n"
                    "#/gate/source/ga68.2/gps/shape                    Sphere                          # Specifying the shape of the voxelised source\n"
                    "#/gate/source/ga68.2/gps/radius                   301.0 um                         # Radius of source\n" 
                    "#/gate/source/ga68.2/gps/radius0                  300.0 um                         # Radius of source\n" 
                    "/gate/source/ga68.2/setForcedUnstableFlag        true                             # Forcing the source to be unstable and have the ability to decay\n"
                    "/gate/source/ga68.2/setForcedHalfLife            4062.6 s                         # Sepecifying the half life of the source, 2.341e+7\n"
                    "/gate/source/ga68.2/gps/centre                   0.0 0.0 0.0 cm                   # Specifying the position of the source\n"
                    "/gate/source/ga68.2/gps/angtype                  iso                              # Specfying the angular distribution of the decay particles\n"
                    "/gate/source/ga68.2/attachTo                     {}                          # Attaching the source to a volume.  This is for motion.  The location os the source is now relative to the volume it is attached to.\n"
                    "\n".format(sourceActivity2, sourceAttach2)
                )

            if f18source2 == True:
                file.write(
                    "#   FLUORINE - 18\n"
                    "/gate/source/addSource                         f18.2 gps                         # ++ Same as Above ++\n"
                    "/gate/source/f18.2/setActivity                   {}                      # ++ Same as Above ++\n" 
                    "/gate/source/f18.2/gps/particle                  e+                              # ++ Same as Above ++\n" 
                    "/gate/source/f18.2/gps/energytype                UserSpectrum                    # ++ Same as Above ++\n" 
                    "/gate/source/f18.2/gps/setSpectrumFile           data/f18_spectrum.txt           # ++ Same as Above ++\n" 
                    "/gate/source/f18.2/setForcedUnstableFlag         true                            # ++ Same as Above ++\n"  
                    "/gate/source/f18.2/setForcedHalfLife             6586.26 s                       # ++ Same as Above ++\n"
                    "/gate/source/f18.2/gps/centre                    0.0 0.0 0.0 cm                  # ++ Same as Above ++\n"
                    "/gate/source/f18.2/gps/angtype                   iso                             # ++ Same as Above ++\n"
                    "/gate/source/f18.2/attachTo                      {}                         # ++ Same as Above ++\n"
                    "\n".format(sourceActivity2, sourceAttach2)
                )

            if na22source2 == True:
                file.write(
                    "#   SODIUM - 22\n"
                    "/gate/source/addSource                          na22.2 gps                        # ++ Same as Above ++\n"
                    "/gate/source/na22.2/setActivity                   {}                      # ++ Same as Above ++\n"
                    "/gate/source/na22.2/gps/particle                  e+                              # ++ Same as Above ++\n"
                    "/gate/source/na22.2/gps/energytype                UserSpectrum                    # ++ Same as Above ++\n"
                    "/gate/source/na22.2/gps/setSpectrumFile           data/na22_spectrum.txt          # ++ Same as Above ++\n"
                    "/gate/source/na22.2/setForcedUnstableFlag         true                            # ++ Same as Above ++\n"
                    "/gate/source/na22.2/setForcedHalfLife             8.199e+7 s                      # ++ Same as Above ++\n"
                    "/gate/source/na22.2/gps/centre                    0.0 0.0 0.0 cm                  # ++ Same as Above ++\n"
                    "/gate/source/na22.2/gps/angtype                   iso                             # ++ Same as Above ++\n"
                    "/gate/source/na22.2/attachTo                      {}                         # ++ Same as Above ++\n"
                    "\n".format(sourceActivity2, sourceAttach2)
                )


        if voxelSource == True:
            file.write(
                "/gate/source/addSource                                             bubbleTest voxel\n"
                "/gate/source/bubbleTest/reader/insert                              image\n"
                "/gate/source/bubbleTest/imageReader/translator/insert              range\n"
                "/gate/source/bubbleTest/imageReader/rangeTranslator/readTable      {}\n"
                "/gate/source/bubbleTest/imageReader/readFile                       {}\n"
                "/gate/source/bubbleTest/imageReader/verbose                        0\n"
                "/gate/source/bubbleTest/setPosition                                {}\n"
                "/gate/source/bubbleTest/dump                                       0\n"
                "/gate/source/bubbleTest/setType                                    {}\n"
                "/gate/source/bubbleTest/gps/particle                               {}\n"
                "/gate/source/bubbleTest/gps/energytype                             {}\n"
                "/gate/source/bubbleTest/gps/monoenergy                             {}\n"
                "/gate/source/bubbleTest/gps/angtype                                {}\n"
                "/gate/source/bubbleTest/gps/confine                                NULL"
                "\n".format(voxelSourcePixel2Activity, voxelSourceImage, voxelSourcePosition, voxelSourceType, voxelSourceParticle, voxelSourceParticleEnergyType, voxelSourceParticleEnergy, voxelSourceAngularDist)
            )

        file.write(
            "/gate/source/list                                                       # Lists the sources being used\n"
            "\n"
        )

        file.write(
            "#   ACQUISITION SETTINGS\n"
            "/gate/application/setTimeSlice     {}                              # Set the time slice, for moving sources controls the granularity of motion\n"
            "/gate/application/setTimeStart     {}   s                               # The acquisition start time\n" 
            "/gate/application/setTimeStop      {}  s                               # The acquisition end time\n"
            "/gate/application/startDAQ                                              # Start the data acquisiion/simulation\n"
            "\n".format(timeStep, coreCounterStart, coreCounterEnd)
        )

        file.close()
        coreCounter += coreTime

    #############################
    ###### MULTIPROCESSING ######
    #############################

    onlyFiles = np.array(os.listdir(MACRO_FOLDER)) 

    def macroRun(path, macro):                                                                                              # Function that will create processes that will be run simualtaneously
        output = run(["Gate", "{}{}".format(path, macro)], stdout = PIPE, universal_newlines = True).stdout
        print(output)

    def atoi(text):
        return int(text) if text.isdigit() else text

    def natural_keys(text):
        return [atoi(c) for c in re.split(r'(\d+)', text)]

    mycwd = os.getcwd()
    os.chdir(MACRO_FOLDER)
    os.system("chmod -R 775 .")                                                                                             # Making it possibe to use the macro files in multiprocessing
    os.chdir(mycwd)
    os.system("chmod -x GATE_multithread_server.py")                                                                        # Making the python script executable, might be unesscary

    process = []                                                                                                            # Store out processes
    for i in range(cores):                                                                                                  # Create the processes
        p = multiprocessing.Process(target = macroRun, args = (MACRO_FOLDER, onlyFiles[i]) )
        p.start()
        process.append(p)

    for j in process:                                                                                                       # Start the processes
        j.join()

    rootFiles = os.listdir(ROOT_FOLDER)                  
    rootFiles.sort(key = natural_keys)                                                                                      # Folder of the root Files
    haddStart = "hadd -f {}".format(OUTPUT_ROOT+OUTPUT_ROOT_FILE_NAME)
    ROOT_FOLDER_LIST = cores * [ROOT_FOLDER]
    a = " ".join([a + b for a,b in zip(ROOT_FOLDER_LIST, rootFiles)])                                                       # String with hadd command
    haddCMD = haddStart + " " + a

    os.system(haddCMD)                                                                                                      # Running the hadd command to combine the root files

    shutil.rmtree(MACRO_FOLDER)                                                                                             # Removes the temporary folders
    shutil.rmtree(ROOT_FOLDER)
