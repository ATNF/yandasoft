# @file
from askapdev.rbuild.builders import Builder
from askapdev.rbuild.utils import run

class MyBuilder(Builder):
    def _build(self):
        pass
        #for fn in ['measdata', 'testdataset']:
        #    run("tar -xjf %s.tar.bz2" % fn)

builder = MyBuilder(pkgname=".")
#builder.add_extra_clean_targets('data')
#builder.add_extra_clean_targets('testdataset.ms')

builder.build()
