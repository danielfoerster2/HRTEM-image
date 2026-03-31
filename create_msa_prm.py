#!/usr/bin/python3

import os
import numpy as np

electron_energy = float(os.getenv('electron_energy'))
nz = int(os.getenv('nz'))

with open('msa.prm', 'w') as prm_file:
    print("'[Microscope Parameters]'", file=prm_file)
    print(f"{25.0}                          (STEM probe forming aperture: radius (mrad), rel. asymmetry, asym. dir. (rad), rel. edge)", file=prm_file)
    print(f"{80.0}                          (lower semi angle of image-space detector (mrad))", file=prm_file)
    print(f"{220.0}                         (upper semi angle of image-space detector (mrad))", file=prm_file)
    print(f"{0} 'detectors.prm'             (switch for using a detector definition file, and the name of the detector definition file, attention: the output file name will change when using more than one detector, the detector definitions in the preceeding two lines are ignored when using a detector definition file, all detector definitions are ignord in CTEM mode)", file=prm_file)
    print(f"{electron_energy}               (electron wavelength (nm), <=1.0, alternatively, >1.0: electron energy (keV))", file=prm_file)
    print(f"{0.0}                           (de-magnified source radius (nm) for applying partial spatial coherence to STEM images)", file=prm_file)
    print(f"{0.0}                           (defocus spread (nm) (1/e half width) for STEM image simulation only)", file=prm_file)
    print(f"{0.0}                           (defocus spread kernel width for explicit focal convolution, STEM only)", file=prm_file)
    print(f"{0}                             (defocus spread kernel steps/size, STEM only)", file=prm_file)
    print(f"{0}                             ([NOA] = number of aberration definitions following this line. Can be zero when no aberration definitions are given.)", file=prm_file)

    print("'[Multislice Parameters]'", file=prm_file)
    print(f"{0.0}                   (object tilt x [deg], approximative approach, keep values smaller than 5 degrees)", file=prm_file)
    print(f"{0.0}                   (object tilt y [deg], approximative approach, keep values smaller than 5 degrees)", file=prm_file)
    print(f"{0.0}                   (horizontal scan frame offset [nm])", file=prm_file)
    print(f"{0.0}                   (vertical scan frame offset [nm])", file=prm_file)
    print(f"{0.0}                   (horizontal scan frame size [nm])", file=prm_file)
    print(f"{0.0}                   (vertical scan frame size [nm])", file=prm_file)
    print(f"{0.0}                   (frame rotation w.r.t. slice x-axis [deg])", file=prm_file)
    print(f"{0}                     (number of scan columns = number of pixels on horizontal scan image axis)", file=prm_file)
    print(f"{0}                     (number of scan rows = number of pixels on vertical scan image axis)", file=prm_file)
    print(f"{0}                     (Switch partial temporal coherence ON = 1 or OFF = 0. Simulation by explicit focal averaging, leadingto a drastic increase of calculation time!)", file=prm_file)
    print(f"{0}                     (Switch partial spatial coherence calculation for an input STEM image. OFF = 0, Gaussian profile = 1, Cauchy profile = 2, Top-hat profile = 3.)", file=prm_file)
    print(f"{1}                     (integer factor to internally repeat supercell data in horizontal direction, x)", file=prm_file)
    print(f"{1}                     (integer factor to internally repeat supercell data in vertical direction, y)", file=prm_file)
    print(f"{1}                     (integer factor to internally repeat supercell data in depth, z, historic obsolete parameter.)", file=prm_file)
    print(f"slice                   (slice file name string [SFN], slice files will be searched with names [SFN]+\"_###.sli\" where ### is a three digit number identifying the indiviual slices in numbered order from entrance slice 001 to exit slice [NOSD]. )", file=prm_file)
    print(f"{nz}                    ([NOSD] = number of slice files in z-order, defines the index range (001 ... [NOSD]) of the slice-file names)", file=prm_file)
    print(f"{1}                     [NFLV] = number of frozen lattice alternatives, if this value is larger than 1. Slice variants can be present either as multiple slice data contained by one file per slice (default), or in form of several individual slice files containing one variation. If the slice variants are present in individual files, the slice file naming is change to [SFN]+\"_###_###.sli\", where the first index identifies the slice variant number from 001 to [NFLV] and the second index identifies the slice number in z-order from 001 to [NOSD].)", file=prm_file)
    print(f"{1}                     ([NCPP] = number of minimum frozen lattice variations for one scan pixel in STEM mode, and number of frozen lattice variations generating exit plane waves in CTEM mode, drastic increase of calculation time when >1, since multiple multi-slice calculations will be done per scan pixel. )", file=prm_file)
    print(f"{0}                     (periodic detector readout positions, number of slices after which detectors are read out. A sequence of output files will be generated, adding the index \"_t###\" where ### is a 3 digit number specifying the slice number of the detection. Default value: 0 -> detectors are read out at the exit plane and just one output file is generated. )", file=prm_file)
    print(f"{nz}                    ([NOS] = total number of slices in the object, minimum this number of slice IDs MUST FOLLOW BELOW! Only the next NOSD lines will be read, additional lines are ignored.)", file=prm_file)
    for i in range(nz):
        print(i, file=prm_file)
