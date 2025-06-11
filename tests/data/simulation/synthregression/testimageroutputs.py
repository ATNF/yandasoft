# regression tests with imager checking outputs produced
# some fixed parameters are given in imageroutputs_template.in

import os.path
from synthprogrunner import *

def checkExists(product, wantIt, spectral):
    postfix = 'iot.taylor.0'
    if spectral:
        postfix = 'iot'

    filename = product+"."+postfix
    if product == "restored":
        if spec:
            filename = "image.restored."+postfix
        else:
            filename = "image."+postfix+".restored"

    if product == "alt.restored":
        filename = "image."+postfix+".alt.restored"
    if product == "alt.residual":
        filename = "residual."+postfix+".alt"

    exists = os.path.exists(filename)

    if wantIt and exists:
        print("INFO found requested output "+filename)

    if not wantIt and not exists:
        print("INFO ouput "+filename+" absent as requested")

    if wantIt and not exists:
        if product == "mask":
            print("INFO Image "+filename+" requested, but not found, this means it was 1.0 everywhere")
        else:
            raise RuntimeError("Image "+filename+" requested, but does not exist")
    if not wantIt and exists:
        raise RuntimeError("Image "+filename+" not requested, but does exist")


def analyseResult(spr, spec, writeResidual, writeRawPsf, writeFirstRestore, restore, writePsfImage,
    writeWtImage, writeMaskImage, writeModelImage, writeSensitivityImage, writeGrids):
   '''
      spr - synthesis program runner (to run imageStats)

      throws exceptions if something is wrong, otherwise just
      returns
   '''

   # Now check for data products
   checkExists("psf", writeRawPsf, spec)
   if (spec):
       checkExists("residual", writeResidual, spec)
       checkExists("restored", restore, True)
   else:
       checkExists("alt.residual", writeResidual, spec)
       checkExists("restored", writeFirstRestore, False)
       checkExists("alt.restored", restore, False)

   checkExists("psf.image", writePsfImage, spec)
   checkExists("weights", writeWtImage, spec)
   checkExists("image", writeModelImage, spec)
   if not spec:
       checkExists("mask", writeMaskImage, False)
       checkExists("sensitivity", writeSensitivityImage, False)
   if spec:
       checkExists("visgrid", writeGrids, True)
       checkExists("pcfgrid", writeGrids, True)
       checkExists("psfgrid", writeGrids, True)

import os

os.system("rm -rf iot.ms")
spr = SynthesisProgramRunner(template_parset = 'simulator.in')
spr.addToParset("Csimulator.dataset = iot.ms")
spr.runSimulator()

for spec in [True, False]:
    for writeIt in [True, False]:
        print("INFO Running imager test : spectral line = "+str(spec)+" write = "+str(writeIt))
        os.system("rm -rf *.iot*")
        writeResidual = writeIt
        writeRawPsf = writeIt
        writeFirstRestore = writeIt
        restore = True
        writePsfImage = writeIt
        writeWtImage = writeIt
        writeMaskImage = writeIt
        writeModelImage = writeIt
        writeSensitivityImage = writeIt
        writeGrids = writeIt

        if spec:
            print("INFO set spectral imaging parameters")
            spr2 = SynthesisProgramRunner(template_parset = 'imager-outputs-spec.in')
        else:
            print("INFO set continuum imaging parameters")
            spr2 = SynthesisProgramRunner(template_parset = 'imager-outputs.in')

        spr2.addToParset("Cimager.dataset = iot.ms")
        spr2.addToParset("Cimager.write.residualimage = "+str(writeResidual))
        spr2.addToParset("Cimager.write.psfrawimage = "+str(writeRawPsf))
        spr2.addToParset("Cimager.write.firstrestore = "+str(writeFirstRestore))
        spr2.addToParset("Cimager.write.psfimage = "+str(writePsfImage))
        spr2.addToParset("Cimager.write.weightsimage = "+str(writeWtImage))
        spr2.addToParset("Cimager.write.maskimage = "+str(writeMaskImage))
        spr2.addToParset("Cimager.write.modelimage = "+str(writeModelImage))
        spr2.addToParset("Cimager.write.sensitivityimage = "+str(writeSensitivityImage))
        spr2.addToParset("Cimager.write.grids = "+str(writeGrids))
        print("INFO About to Run new Imager")

        if "CI" in os.environ:
            spr2.runNewImagerParallel(4)
        else:
            spr2.runNewImagerParallel(4)

        analyseResult(spr2, spec, writeResidual, writeRawPsf, writeFirstRestore, restore, writePsfImage,
            writeWtImage, writeMaskImage, writeModelImage, writeSensitivityImage, writeGrids)

# clean up
os.system("rm -rf iot.ms *.iot* temp_parset.in")
