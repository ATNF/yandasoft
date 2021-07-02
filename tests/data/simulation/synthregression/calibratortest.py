# regression tests of calibrator
# some fixed parameters are given in calibratortest_template.in

import argparse
import cmath

from synthprogrunner import *

#from askapdev.rbuild import setup
#from askapdev.rbuild.dependencies import Dependency

#dep = Dependency(silent=False)
#dep.DEPFILE = "../../dependencies"
#dep.add_package()

#import askap.parset

def analyseResult(spr):
   '''
      spr - synthesis program runner (to run imageStats)

      throws exceptions if something is wrong, otherwise just
      returns
   '''
   src_offset = 0.004/math.pi*180.
   psf_peak=[-172.5,-45]
   true_peak=sinProjection(psf_peak,src_offset,src_offset)
   stats = spr.imageStats('image.field1.restored')
   print("Statistics for restored image: ",stats)
   disterr = getDistance(stats,true_peak[0],true_peak[1])*3600.
   if disterr > 8:
      raise RuntimeError("Offset between true and expected position exceeds 1 cell size (8 arcsec), d=%f, true_peak=%s" % (disterr,true_peak))
   # as polarisation conversion is now fixed in the component-based measurement equation we have exactly the same flux value as in the simulation parset
   if abs(stats['peak']-1.)>0.1:
      raise RuntimeError("Peak flux in the image is notably different from 1 Jy, F=%f" % stats['peak'])

   stats = spr.imageStats('image.field1')
   print("Statistics for modelimage: ",stats)
   disterr = getDistance(stats,true_peak[0],true_peak[1])*3600.
   if disterr > 8:
      raise RuntimeError("Offset between true and expected position exceeds 1 cell size (8 arcsec), d=%f, true_peak=%s" % (disterr,true_peak))

   stats = spr.imageStats('psf.field1')
   print("Statistics for psf image: ",stats)
   disterr = getDistance(stats,psf_peak[0],psf_peak[1])*3600.
   if disterr > 8:
      raise RuntimeError("Offset between true and expected position exceeds 1 cell size (8 arcsec), d=%f, true_peak=%s" % (disterr,true_peak))

   stats = spr.imageStats('weights.field1')
   print("Statistics for weight image: ",stats)
   if abs(stats['rms']-stats['peak'])>0.1 or abs(stats['rms']-stats['median'])>0.1 or abs(stats['peak']-stats['median'])>0.1:
      raise RuntimeError("Weight image is expected to be constant for WProject and WStack gridders")

   stats = spr.imageStats('residual.field1')
   print("Statistics for residual image: ",stats)
   if stats['rms']>0.01 or abs(stats['median'])>0.0001:
      raise RuntimeError("Residual image has too high rms or median. Please verify")


def loadParset(fname,rotate=True):
    """
       Helper method to read parset file with gains into a python dictionary. Most likely we wouldn't need this method
       if importing of py-parset was a bit more straightforward. 

       fname - file name of the parset file to read
       rotate - if True, phases of all gains will be rotated, so the first antenna will get 0 phase on both polarisations

       Return: dictionary with gains (name is the key, value is the complex value).

       For simplicity we assume that both real and imaginary parts are given as this is the case for all parset files
       this helper method is likely to be used for.
    """
    res = {}
    f = open(fname)
    try:
       for line in f:
          if line.startswith("gain"):
             parts = line.split()
             if len(parts)!=3:
                 raise RuntimeError("Expect 3 parts in the line of parset file, you have %s" % (parts,))
             if parts[1]!="=":
                 raise RuntimeError("Value and key are supposed to be separated by '=', you have %s" % (parts,))
             if not parts[2].startswith('[') or not parts[2].endswith(']'):
                 raise RuntimeErrror("Value is supposed to be in square brackets, you have <%s>" % parts[2])
             values = parts[2][1:-1].split(",")
             if len(values)!=2:
                 raise RuntimeError("Two numbers are expected, you have %s" % (values,))
             res[parts[0]] = float(values[0])+(1j)*float(values[1])
    finally:
       f.close()
    if rotate:
       first_pol = 1.
       second_pol = 1.
       if "gain.g11.0.0" in res:
          first_pol = res["gain.g11.0.0"].conjugate()
       first_pol /= abs(first_pol)
       if "gain.g22.0.0" in res:
          second_pol = res["gain.g22.0.0"].conjugate()
       second_pol /= abs(first_pol)
       for k,v in list(res.items()):
          if "gain.g11" in k:
              res[k] *= first_pol
          if "gain.g22" in k:
              res[k] *= second_pol

    return res

def runTests(solverType):
    """
    Runs tests for a given calibration solver type.
    """
    spr = SynthesisProgramRunner(template_parset = 'calibratortest_template.in')
    spr.addToParset("Csimulator.corrupt = false")
    spr.runSimulator()
    
    spr.initParset()
    spr.runImager()
    analyseResult(spr)
    
    print("First run of ccalibrator, should get gains close to (1.,0.)")
    
    spr.addToParset("Ccalibrator.calibaccess = parset")
    spr.addToParset("Ccalibrator.solver = " + solverType)
    spr.runCalibrator()
    # here result.dat should be close to (1.,0) within 0.03 or so
    
    res_gains = loadParset("result.dat")
    for k,v in list(res_gains.items()):
       if abs(v-1)>0.03:
          raise RuntimeError("Gain parameter %s has a value of %s which is notably different from (1,0)" % (k,v))
    
    # now repeat the simulation, but with corruption of visibilities
    spr.initParset()
    spr.addToParset("Csimulator.corrupt = true")
    spr.runSimulator()
    
    print("Second run of ccalibrator, gains should be close to rndgains.in")
    # calibrate again
    spr.addToParset("Ccalibrator.calibaccess = parset")
    spr.addToParset("Ccalibrator.solver = " + solverType)
    spr.runCalibrator()
    
    # gains should now be close to rndgains.in
    
    res_gains = loadParset("result.dat")
    orig_gains = loadParset("rndgains.in")
    for k,v in list(res_gains.items()):
       if k not in orig_gains:
          raise RintimeError("Gain parameter %s found in the result is missing in the model!" % k)
       orig_val = orig_gains[k]
       if abs(v-orig_val)>0.03:
          raise RuntimeError("Gain parameter %s has a value of %s which is notably different from model value %s" % (k,v,orig_val))
    
    # now try to obtain time-dependent solution (note, a proper analysis of the result is not done)
    print("Third run of ccalibrator. Time-dependent solution for antennagains")
    spr.initParset()
    spr.addToParset("Ccalibrator.calibaccess = table")
    spr.addToParset("Ccalibrator.interval = 600s")
    spr.addToParset("Ccalibrator.solve = antennagains")
    spr.addToParset("Ccalibrator.solver = " + solverType)
    os.system("rm -rf caldata.tab")
    spr.runCalibrator()
    
    print("Testing calibration application.")
    # run cimager applying time-dependent calibration
    spr.addToParset("Cimager.calibrate = true")
    spr.addToParset("Cimager.calibaccess = table")
    spr.addToParset("Cimager.calibaccess.table = \"caldata.tab\"")
    spr.addToParset("Cimager.calibrate.ignorebeam = true")
    spr.runImager()
    analyseResult(spr)

def runTestsParallel(ncycles):
    """
    Runs parallel tests for bandpass calibration using ccalibrator app, with parallel matrix using LSQR solver.
    The test is performed by running serial and parallel versions and comparing the results.
    Basically, testing that a parallel matrix case (one big matrix with all channels) produces
    the same results as serial matrices case (several small matrixes for each channel chunk).
    """
    msarchive = "vis_4chan.tar.bz2"
    msfile = "vis_4chan.ms"

    # Extracting measurement set from an archive.
    if not os.path.exists(msarchive):
        raise RuntimeError("A tarball with measurement sets does not seem to exist (%s)" % msarchive)

    if os.path.exists(msfile):
        print("Removing old %s" % msfile)
        os.system("rm -rf %s" % msfile)

    os.system("tar -xjf %s" % msarchive)

    #----------------------------------------------------------------------------------
    tmp_parset = 'ccal_tmp.in'
    os.system("rm -f %s" % tmp_parset)
    os.system("touch %s" % tmp_parset)

    spr = SynthesisProgramRunner(template_parset = tmp_parset)

    spr.addToParset("Ccalibrator.dataset                      = [%s]" % msfile)
    spr.addToParset("Ccalibrator.nAnt                         = 36")
    spr.addToParset("Ccalibrator.nBeam                        = 1")
    spr.addToParset("Ccalibrator.nChan                        = 4")

    spr.addToParset("Ccalibrator.calibaccess                  = parset")
    spr.addToParset("Ccalibrator.calibaccess.parset           = result.dat")

    spr.addToParset("Ccalibrator.sources.names                = [field1]")
    spr.addToParset("Ccalibrator.sources.field1.direction     = [12h30m00.000, -45.00.00.000, J2000]")
    spr.addToParset("Ccalibrator.sources.field1.components    = [src1]")
    spr.addToParset("Ccalibrator.sources.src1.flux.i          = 1.0")
    spr.addToParset("Ccalibrator.sources.src1.direction.ra    = 0.00")
    spr.addToParset("Ccalibrator.sources.src1.direction.dec   = 0.00")

    spr.addToParset("Ccalibrator.gridder                      = SphFunc")

    spr.addToParset("Ccalibrator.ncycles                      = %s" % ncycles)

    spr.addToParset("Ccalibrator.solve                        = bandpass")
    spr.addToParset("Ccalibrator.solver                       = LSQR")
    spr.addToParset("Ccalibrator.solver.LSQR.verbose          = true")

    # Filenames for writing the calibration results (parsets).
    result_serial = "result_serial.dat"
    result_parallel = "result_parallel.dat"

    # Serial run.
    print("Bandpass test: Serial run of ccalibrator.")
    nprocs = 1
    spr.runCalibratorParallel(nprocs)

    # Store the results.
    os.system("mv result.dat %s" % result_serial)

    # Set additional parameters for the parallel run.
    # These are defined for the case of 4 channels with 5 cpus, i.e., 1 channel per worker (exluding master rank).
    spr.addToParset("Ccalibrator.chanperworker                = 1")
    spr.addToParset("Ccalibrator.chunk                        = %w")
    spr.addToParset("Ccalibrator.solver.LSQR.parallelMatrix   = true")

    # Parallel run.
    print("Bandpass test: Parallel run of ccalibrator.")
    nprocs = 5
    spr.runCalibratorParallel(nprocs)

    # Store the results.
    os.system("mv result.dat %s" % result_parallel)

    tol = 1.e-6
    compareGains(result_serial, result_parallel, tol)

def runTestsSmoothnessConstraintsParallel():
    """
    Runs parallel tests for bandpass calibration using ccalibrator app, using LSQR solver with smoothness constraints.
    The test is performed by running on one and several workers, and then comparing the output gains (expected to be the same).
    """
    msarchive = "40chan.ms.tar.bz2"
    msfile = "40chan.ms"

    # Extracting measurement set from an archive.
    if not os.path.exists(msarchive):
        raise RuntimeError("A tarball with measurement sets does not seem to exist (%s)" % msarchive)

    if os.path.exists(msfile):
        print("Removing old %s" % msfile)
        os.system("rm -rf %s" % msfile)

    os.system("tar -xjf %s" % msarchive)

    #----------------------------------------------------------------------------------
    tmp_parset = 'ccal_tmp.in'
    os.system("rm -f %s" % tmp_parset)
    os.system("touch %s" % tmp_parset)

    spr = SynthesisProgramRunner(template_parset = tmp_parset)

    spr.addToParset("Ccalibrator.dataset                      = [%s]" % msfile)
    spr.addToParset("Ccalibrator.nAnt                         = 12")
    spr.addToParset("Ccalibrator.nBeam                        = 1")
    spr.addToParset("Ccalibrator.nChan                        = 40")
    
    spr.addToParset("Ccalibrator.refantenna                   = 1")

    spr.addToParset("Ccalibrator.calibaccess                  = parset")
    spr.addToParset("Ccalibrator.calibaccess.parset           = result.dat")

    spr.addToParset("Ccalibrator.sources.names                = [field1]")
    spr.addToParset("Ccalibrator.sources.field1.direction     = [19h39m25.036, -63.42.45.63, J2000]")
    spr.addToParset("Ccalibrator.sources.field1.components    = [src]")
    spr.addToParset("Ccalibrator.sources.src.calibrator       = 1934-638")

    spr.addToParset("Ccalibrator.gridder                      = SphFunc")

    spr.addToParset("Ccalibrator.ncycles                      = 20")

    spr.addToParset("Ccalibrator.solve                        = bandpass")
    spr.addToParset("Ccalibrator.solver                       = LSQR")
    spr.addToParset("Ccalibrator.solver.LSQR.verbose          = true")
    spr.addToParset("Ccalibrator.solver.LSQR.smoothing        = true")
    spr.addToParset("Ccalibrator.solver.LSQR.alpha            = 1e5")
    spr.addToParset("Ccalibrator.solver.LSQR.parallelMatrix   = true")

    spr.addToParset("Ccalibrator.chunk                        = %w")

    # Filenames for writing the calibration results (parsets).
    result_w1 = "result_smoothing_w1.dat"
    result_w4 = "result_smoothing_w4.dat"
    result_w10 = "result_smoothing_w10.dat"

    #------------------------------------------------------------------------------
    # Set partitioning for 40 channels with 2 cpus, i.e., all 40 channels per worker (exluding master rank).
    spr.addToParset("Ccalibrator.chanperworker                = 40")

    print("Bandpass test: One worker.")
    nprocs = 2
    spr.runCalibratorParallel(nprocs)

    # Store the results.
    os.system("mv result.dat %s" % result_w1)

    #------------------------------------------------------------------------------
    # Set partitioning for 40 channels with 5 cpus, i.e., 10 channels per worker (exluding master rank).
    spr.addToParset("Ccalibrator.chanperworker                = 10")

    print("Bandpass test: Four workers.")
    nprocs = 5
    spr.runCalibratorParallel(nprocs)

    # Store the results.
    os.system("mv result.dat %s" % result_w4)

    #------------------------------------------------------------------------------
    # Set partitioning for 40 channels with 11 cpus, i.e., 4 channels per worker (exluding master rank).
    spr.addToParset("Ccalibrator.chanperworker                = 4")

    # Parallel run.
    print("Bandpass test: Ten workers.")
    nprocs = 11
    spr.runCalibratorParallel(nprocs)

    # Store the results.
    os.system("mv result.dat %s" % result_w10)

    #------------------------------------------------------------------------------
    # Compare the output results.
    tol = 1.e-4
    compareGains(result_w1, result_w4, tol)
    compareGains(result_w1, result_w10, tol)
    compareGains(result_w4, result_w10, tol)

def runTestsSmoothnessConstraintsGradientCost():
    """
    Runs tests for bandpass calibration using ccalibrator app, using LSQR solver with smoothness constraints.
    The test is performed by running calibration with and without smoothness constraints and comparing the gradient cost to expected values.
    """
    msarchive = "40chan.ms.tar.bz2"
    msfile = "40chan.ms"

    # Extracting measurement set from an archive.
    if not os.path.exists(msarchive):
        raise RuntimeError("A tarball with measurement sets does not seem to exist (%s)" % msarchive)

    if os.path.exists(msfile):
        print("Removing old %s" % msfile)
        os.system("rm -rf %s" % msfile)

    os.system("tar -xjf %s" % msarchive)

    #----------------------------------------------------------------------------------
    tmp_parset = 'ccal_tmp.in'
    os.system("rm -f %s" % tmp_parset)
    os.system("touch %s" % tmp_parset)

    spr = SynthesisProgramRunner(template_parset = tmp_parset)

    spr.addToParset("Ccalibrator.dataset                      = [%s]" % msfile)
    spr.addToParset("Ccalibrator.nAnt                         = 12")
    spr.addToParset("Ccalibrator.nBeam                        = 1")
    spr.addToParset("Ccalibrator.nChan                        = 40")
    
    spr.addToParset("Ccalibrator.refantenna                   = 1")

    spr.addToParset("Ccalibrator.calibaccess                  = parset")
    spr.addToParset("Ccalibrator.calibaccess.parset           = result.dat")

    spr.addToParset("Ccalibrator.sources.names                = [field1]")
    spr.addToParset("Ccalibrator.sources.field1.direction     = [19h39m25.036, -63.42.45.63, J2000]")
    spr.addToParset("Ccalibrator.sources.field1.components    = [src]")
    spr.addToParset("Ccalibrator.sources.src.calibrator       = 1934-638")

    spr.addToParset("Ccalibrator.gridder                      = SphFunc")

    spr.addToParset("Ccalibrator.ncycles                      = 20")

    spr.addToParset("Ccalibrator.solve                        = bandpass")
    spr.addToParset("Ccalibrator.solver                       = LSQR")
    spr.addToParset("Ccalibrator.solver.LSQR.verbose          = true")
    spr.addToParset("Ccalibrator.solver.LSQR.smoothing        = false")
    spr.addToParset("Ccalibrator.solver.LSQR.alpha            = 1e5")
    spr.addToParset("Ccalibrator.solver.LSQR.parallelMatrix   = true")

    spr.addToParset("Ccalibrator.chunk                        = %w")
    spr.addToParset("Ccalibrator.chanperworker                = 40")

    # Set partitioning for 40 channels with 2 cpus, i.e., all 40 channels on one worker (plus a master rank).
    nprocs = 2

    # Filenames for writing the calibration results (parsets).
    result_nonsmooth = "result_nonsmooth.dat"
    result_smooth = "result_smooth.dat"
    result_smooth2 = "result_smooth2.dat"
    result_smooth4 = "result_smooth4.dat"
    result_smooth2_acc4 = "result_smooth2_acc4.dat"
    result_smooth4_acc4 = "result_smooth4_acc4.dat"

    # Remove old results (if any).
    os.system("rm %s" % result_nonsmooth)
    os.system("rm %s" % result_smooth)
    os.system("rm %s" % result_smooth2)
    os.system("rm %s" % result_smooth4)
    os.system("rm %s" % result_smooth2_acc4)
    os.system("rm %s" % result_smooth4_acc4)

    #------------------------------------------------------------------------------
    print("Bandpass test: without smoothing constraints.")
    spr.runCalibratorParallel(nprocs)

    # Store the results.
    os.system("mv result.dat %s" % result_nonsmooth)

    #------------------------------------------------------------------------------
    # Switch on the smoothing constraints (with Gradient smoother).
    spr.addToParset("Ccalibrator.solver.LSQR.smoothing            = true")
    spr.addToParset("Ccalibrator.solver.LSQR.smoothing.maxWeight  = 3.e+6")
    spr.addToParset("Ccalibrator.solver.LSQR.smoothing.type       = 0")

    print("Bandpass test: with smoothing constraints (type = 0).")
    spr.runCalibratorParallel(nprocs)

    # Store the results.
    os.system("mv result.dat %s" % result_smooth)

    #------------------------------------------------------------------------------
    # Switch on the smoothing constraints (with Laplacian smoother).
    spr.addToParset("Ccalibrator.solver.LSQR.smoothing.type       = 2")
    spr.addToParset("Ccalibrator.solver.LSQR.smoothing.accuracy   = 2")

    print("Bandpass test: with smoothing constraints (type = 2).")
    spr.runCalibratorParallel(nprocs)

    # Store the results.
    os.system("mv result.dat %s" % result_smooth2)

    #-------------------------------
    # Use the 4th order of accuracy.
    spr.addToParset("Ccalibrator.solver.LSQR.smoothing.accuracy    = 4")

    print("Bandpass test: with smoothing constraints (type = 2), accuracy 4.")
    spr.runCalibratorParallel(nprocs)

    # Store the results.
    os.system("mv result.dat %s" % result_smooth2_acc4)

    #------------------------------------------------------------------------------
    # Switch on the smoothing constraints (with 4th order smoother).
    spr.addToParset("Ccalibrator.solver.LSQR.smoothing.type       = 4")
    spr.addToParset("Ccalibrator.solver.LSQR.smoothing.accuracy   = 2")

    print("Bandpass test: with smoothing constraints (type = 4).")
    spr.runCalibratorParallel(nprocs)

    # Store the results.
    os.system("mv result.dat %s" % result_smooth4)

    #-------------------------------
    # Use the 4th order of accuracy.
    spr.addToParset("Ccalibrator.solver.LSQR.smoothing.accuracy    = 4")

    print("Bandpass test: with smoothing constraints (type = 4), accuracy 4.")
    spr.runCalibratorParallel(nprocs)

    # Store the results.
    os.system("mv result.dat %s" % result_smooth4_acc4)

    #------------------------------------------------------------------------------
    # Calcualte the gradient cost.
    nchan = 40
    nant = 12
    cost_nonsmooth = calculateGradientCost(result_nonsmooth, nchan, nant, True)
    cost_smooth = calculateGradientCost(result_smooth, nchan, nant, True)

    # Note, for Laplacian/4th-order cases we also calculate gradient here (instead of Laplacian/4th-order),
    # just to compare results using a single number, so that the test fails if we break the code.
    cost_smooth2 = calculateGradientCost(result_smooth2, nchan, nant, True)
    cost_smooth4 = calculateGradientCost(result_smooth4, nchan, nant, True)
    cost_smooth2_acc4 = calculateGradientCost(result_smooth2_acc4, nchan, nant, True)
    cost_smooth4_acc4 = calculateGradientCost(result_smooth4_acc4, nchan, nant, True)

    print('cost nonsmooth =', cost_nonsmooth)
    print('cost smooth =', cost_smooth)
    print('cost smooth2 =', cost_smooth2)
    print('cost smooth4 =', cost_smooth4)
    print('cost smooth2 acc4 =', cost_smooth2_acc4)
    print('cost smooth4 acc4 =', cost_smooth4_acc4)

    #------------------------------------------------------------------------------
    # Verify the gradient cost.

    # Note that the expected_cost_smooth value calculated here differs from the one printed in the log (~0.0545)
    # due to phase referencing performed in the end of calibration, which does not preserve this type of the gradient.
    expected_cost_nonsmooth = 22.645816
    expected_cost_smooth = 0.078488
    expected_cost_smooth2 = 0.188401
    expected_cost_smooth4 = 0.666211
    expected_cost_smooth2_acc4 = 0.833353
    expected_cost_smooth4_acc4 = 1.393525

    tol = 1.e-3
    if abs(cost_nonsmooth - expected_cost_nonsmooth) > tol:
        raise RuntimeError("Nonsmooth gradient cost is wrong! cost = %s" % cost_nonsmooth)

    tol = 1.e-5
    if abs(cost_smooth - expected_cost_smooth) > tol:
        raise RuntimeError("Smooth gradient cost is wrong! cost = %s" % cost_smooth)

    tol = 1.e-5
    if abs(cost_smooth2 - expected_cost_smooth2) > tol:
        raise RuntimeError("Smooth gradient2 cost is wrong! cost = %s" % cost_smooth2)

    tol = 1.e-5
    if abs(cost_smooth4 - expected_cost_smooth4) > tol:
        raise RuntimeError("Smooth gradient4 cost is wrong! cost = %s" % cost_smooth4)

    tol = 1.e-5
    if abs(cost_smooth2_acc4 - expected_cost_smooth2_acc4) > tol:
        raise RuntimeError("Smooth gradient2 acc4 cost is wrong! cost = %s" % cost_smooth2_acc4)

    tol = 1.e-5
    if abs(cost_smooth4_acc4 - expected_cost_smooth4_acc4) > tol:
        raise RuntimeError("Smooth gradient4 acc4 cost is wrong! cost = %s" % cost_smooth4_acc4)

def compareGains(file1, file2, tol):
    gains1 = loadParset(file1)
    gains2 = loadParset(file2)

    # Comparing the results between the serial and parallel runs.
    for k, v in list(gains1.items()):
        if k not in gains2:
            raise RintimeError("Gain parameter %s is missing!" % k)
        val1 = gains1[k]
        val2 = gains2[k]

        if abs(val1 - val2) > tol:
            raise RuntimeError("Gain parameter %s has a value value of %s which is different from value %s" % (k, val1, val2))

def calculateGradientCost(filename, nchan, nant, useComplexNumberParts):
    gains = loadParset(filename)

    pols = ["g11", "g22"]

    grad_cost = 0.
    for pol in pols:
        for ant in range(0, nant):
            for chan in range(0, nchan - 1):
                base_name = "gain." + pol + "." + str(ant) + ".0."
                curr_parname = base_name + str(chan)        # current channel
                next_parname = base_name + str(chan + 1)    # next channel

                curr_gain_complex_val = gains[curr_parname]
                next_gain_complex_val = gains[next_parname]

                # Calculate gradient in frequency, using forward difference approximation: x[i]' = x[i+1] - x[i]
                if useComplexNumberParts:
                    # Apply gradient directly on X and Y parts of a complex number g = X + iY.
                    grad_cost += (next_gain_complex_val.real - curr_gain_complex_val.real)**2
                    grad_cost += (next_gain_complex_val.imag - curr_gain_complex_val.imag)**2
                else:
                    # Apply gradient on the magnitude and phase of a complex number.
                    curr_gain_magnitude = abs(curr_gain_complex_val)
                    next_gain_magnitude = abs(next_gain_complex_val)

                    curr_gain_phase = cmath.phase(curr_gain_complex_val)
                    next_gain_phase = cmath.phase(next_gain_complex_val)

                    grad_cost += (next_gain_magnitude - curr_gain_magnitude)**2
                    grad_cost += (next_gain_phase - curr_gain_phase)**2

    return grad_cost

#==================================================================================================
if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--parallel', help="Include parallel tests", action='store_true')
    opt = parser.parse_args()

    runTests("SVD")
    runTests("LSQR")
    if opt.parallel:
        runTestsParallel(1)
        runTestsParallel(5)
        runTestsSmoothnessConstraintsParallel()
        runTestsSmoothnessConstraintsGradientCost()
