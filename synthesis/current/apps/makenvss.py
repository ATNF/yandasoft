import math
import pymeasures
dm=pymeasures.measures()

import pycasatable

def putline(t, name, value):
    t.addrows(1);
    row=t.nrows()-1
    t.putcell('NAME', row, name)
    t.putcell('VALUES', row, value)
    t.putcell('FREE', row, False)
    t.putcell('AXES', row, '')
    t.putcell('AXESSTART', row, 0.0)
    t.putcell('AXESEND', row, 0.0)
    t.putcell('DOMAIN', row, '')
    t.putcell('DOMAINSTART', row, 0.0)
    t.putcell('DOMAINEND', row, 0.0)
    
def addsource(t, dir, rad, flux, bmaj, bmin, bpa):
    arcsec=math.pi/(180*3600.0)
    degree=math.pi/180.0
    source='source%04d'%(t.nrows()/6)
    ra=dir['m0']['value']
    dec=dir['m1']['value']
    putline(t, 'flux.i.'+source, flux)
    putline(t, 'direction.ra.'+source, ra)
    putline(t, 'direction.dec.'+source, dec)
    putline(t, 'shape.bmaj.'+source, bmaj*arcsec)
    putline(t, 'shape.bmin.'+source, bmin*arcsec)
    putline(t, 'shape.bpa.'+source, bpa*degree)

tabledesc={'AXES': {'_c_order': True,
                    'comment': '',
                    'dataManagerGroup': 'StandardStMan',
                    'dataManagerType': 'StandardStMan',
                    'maxlen': 0,
                    'ndim': -1,
                    'option': 0,
                    'valueType': 'string'},
           'AXESEND': {'_c_order': True,
                       'comment': '',
                       'dataManagerGroup': 'StandardStMan',
                       'dataManagerType': 'StandardStMan',
                       'maxlen': 0,
                       'ndim': 1,
                       'option': 0,
                       'valueType': 'double'},
           'AXESSTART': {'_c_order': True,
                         'comment': '',
                         'dataManagerGroup': 'StandardStMan',
                         'dataManagerType': 'StandardStMan',
                         'maxlen': 0,
                         'ndim': 1,
                         'option': 0,
                         'valueType': 'double'},
           'DOMAIN': {'_c_order': True,
                      'comment': '',
                      'dataManagerGroup': 'StandardStMan',
                      'dataManagerType': 'StandardStMan',
                      'maxlen': 0,
                      'ndim': -1,
                      'option': 0,
                      'valueType': 'string'},
           'DOMAINEND': {'_c_order': True,
                         'comment': '',
                         'dataManagerGroup': 'StandardStMan',
                         'dataManagerType': 'StandardStMan',
                         'maxlen': 0,
                         'ndim': 1,
                         'option': 0,
                         'valueType': 'double'},
           'DOMAINSTART': {'_c_order': True,
                           'comment': '',
                           'dataManagerGroup': 'StandardStMan',
                           'dataManagerType': 'StandardStMan',
                           'maxlen': 0,
                           'ndim': 1,
                           'option': 0,
                           'valueType': 'double'},
           'FREE': {'comment': '',
                    'dataManagerGroup': 'StandardStMan',
                    'dataManagerType': 'StandardStMan',
                    'maxlen': 0,
                    'option': 0,
                    'valueType': 'boolean'},
           'NAME': {'comment': '',
                    'dataManagerGroup': 'StandardStMan',
                    'dataManagerType': 'StandardStMan',
                    'maxlen': 0,
                    'option': 0,
                    'valueType': 'string'},
           'VALUES': {'_c_order': True,
                      'comment': '',
                      'dataManagerGroup': 'StandardStMan',
                      'dataManagerType': 'StandardStMan',
                      'maxlen': 0,
                      'ndim': -1,
                      'option': 0,
                      'valueType': 'double'},
           '_define_hypercolumn_': {}}

t=pycasatable.table('nvss.par', tabledesc=tabledesc, readonly=False)

addsource(t, dir=dm.direction('J2000', '12h01m53.09', '-44.39.39.2'),
          rad=17958,
          flux=  235.7,
          bmaj=  0.0,
          bmin=  0.0,
          bpa= -46.5)

addsource(t, dir=dm.direction('J2000', '12h02m08.43', '-44.44.20.8'),
          rad=17767,
          flux=  105.0,
          bmaj=  33.5,
          bmin=  17.2,
          bpa=  85.1)

addsource(t, dir=dm.direction('J2000', '12h02m19.54', '-45.13.23.2'),
          rad=17569,
          flux=   13.3,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h02m20.02', '-45.14.45.6'),
          rad=17564,
          flux=   16.6,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h02m59.57', '-44.05.48.4'),
          rad=17598,
          flux=   36.2,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h03m09.28', '-44.31.54.0'),
          rad=17212,
          flux=   94.4,
          bmaj= 107.3,
          bmin=  74.0,
          bpa= -48.7)

addsource(t, dir=dm.direction('J2000', '12h03m26.49', '-44.36.36.5'),
          rad=16993,
          flux=   21.6,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h03m30.45', '-46.55.51.0'),
          rad=17942,
          flux=   12.7,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h03m31.09', '-45.28.54.9'),
          rad=16848,
          flux=   17.5,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h03m31.54', '-45.29.47.3'),
          rad=16846,
          flux=   17.5,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h03m31.75', '-46.32.55.4'),
          rad=17500,
          flux=  416.9,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h03m35.49', '-45.10.49.8'),
          rad=16769,
          flux=   38.5,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h03m36.36', '-43.33.31.5'),
          rad=17753,
          flux=  260.5,
          bmaj=  14.4,
          bmin=  0.0,
          bpa=  50.4)

addsource(t, dir=dm.direction('J2000', '12h03m38.28', '-43.49.06.9'),
          rad=17448,
          flux=   76.0,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h03m44.41', '-46.08.42.8'),
          rad=17025,
          flux=   28.0,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h03m50.21', '-43.45.42.7'),
          rad=17383,
          flux=   15.5,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h03m59.28', '-44.30.09.3'),
          rad=16699,
          flux=   16.9,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h04m01.19', '-44.36.10.1'),
          rad=16630,
          flux=   29.7,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h04m14.07', '-46.16.29.0'),
          rad=16826,
          flux=   51.3,
          bmaj=  21.4,
          bmin=  0.0,
          bpa=   9.4)

addsource(t, dir=dm.direction('J2000', '12h04m19.34', '-45.07.56.6'),
          rad=16308,
          flux=   68.3,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h04m27.63', '-46.52.49.2'),
          rad=17331,
          flux=  141.6,
          bmaj=  28.6,
          bmin=  0.0,
          bpa=  34.1)

addsource(t, dir=dm.direction('J2000', '12h04m39.96', '-47.22.08.7'),
          rad=17912,
          flux=   57.5,
          bmaj=  25.2,
          bmin= 0.0,
          bpa= -82.1)

addsource(t, dir=dm.direction('J2000', '12h04m54.60', '-44.29.54.7'),
          rad=16118,
          flux=   12.3,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h04m59.23', '-45.59.36.2'),
          rad=16158,
          flux=   96.1,
          bmaj=  60.3,
          bmin=  0.0,
          bpa= -18.7)

addsource(t, dir=dm.direction('J2000', '12h05m01.51', '-47.04.02.2'),
          rad=17261,
          flux=   15.8,
          bmaj=  22.8,
          bmin= 0.0,
          bpa= -47.4)

addsource(t, dir=dm.direction('J2000', '12h05m03.11', '-42.59.58.8'),
          rad=17656,
          flux=   71.1,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h05m10.72', '-42.49.10.2'),
          rad=17877,
          flux=   12.7,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h05m14.79', '-46.18.09.5'),
          rad=16242,
          flux=   17.5,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h05m15.64', '-43.10.08.6'),
          rad=17272,
          flux=   21.5,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h05m16.46', '-45.20.37.7'),
          rad=15718,
          flux=   13.3,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h05m23.11', '-43.22.33.6'),
          rad=16903,
          flux=   10.2,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h05m36.95', '-45.03.32.9'),
          rad=15493,
          flux=   40.3,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h05m38.77', '-43.25.52.4'),
          rad=16670,
          flux=   15.9,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h05m44.51', '-44.43.40.5'),
          rad=15487,
          flux=   25.8,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h05m45.38', '-44.03.25.8'),
          rad=15901,
          flux=   17.7,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h05m51.30', '-44.29.10.2'),
          rad=15527,
          flux=   56.4,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h05m52.22', '-47.45.58.1'),
          rad=17958,
          flux=   16.9,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h05m52.38', '-43.42.36.8'),
          rad=16184,
          flux=   36.8,
          bmaj=  30.2,
          bmin=  0.0,
          bpa=  31.2)

addsource(t, dir=dm.direction('J2000', '12h06m01.34', '-47.28.46.0'),
          rad=17363,
          flux=   51.4,
          bmaj= 126.7,
          bmin=  55.3,
          bpa=  16.6)

addsource(t, dir=dm.direction('J2000', '12h06m05.47', '-44.03.18.8'),
          rad=15694,
          flux=   14.5,
          bmaj=  23.2,
          bmin= 0.0,
          bpa= -76.0)

addsource(t, dir=dm.direction('J2000', '12h06m08.31', '-42.42.47.7'),
          rad=17510,
          flux=   12.6,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h06m09.25', '-46.45.31.4'),
          rad=16205,
          flux=   10.2,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h06m09.41', '-45.26.31.8'),
          rad=15181,
          flux=   12.2,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h06m09.42', '-42.30.42.3'),
          rad=17872,
          flux=   18.6,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h06m10.25', '-42.39.02.9'),
          rad=17605,
          flux=   14.0,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h06m15.06', '-45.34.31.2'),
          rad=15162,
          flux=   43.0,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h06m15.12', '-42.30.26.7'),
          rad=17828,
          flux=   65.0,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h06m23.63', '-42.34.38.8'),
          rad=17615,
          flux=   38.9,
          bmaj=  79.4,
          bmin=  62.0,
          bpa= -86.3)

addsource(t, dir=dm.direction('J2000', '12h06m24.43', '-44.54.52.8'),
          rad=15012,
          flux=   24.1,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h06m27.06', '-43.40.46.2'),
          rad=15866,
          flux=   20.3,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h06m31.43', '-45.33.34.1'),
          rad=14986,
          flux=   39.0,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h06m35.95', '-42.41.35.9'),
          rad=17284,
          flux=   16.9,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h06m41.50', '-43.08.10.0'),
          rad=16476,
          flux=   16.6,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h06m45.02', '-46.43.37.1'),
          rad=15822,
          flux=   36.6,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h06m45.62', '-46.10.15.2'),
          rad=15215,
          flux=   17.0,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h06m46.84', '-47.07.40.4'),
          rad=16376,
          flux=   21.6,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h06m53.63', '-44.25.27.7'),
          rad=14907,
          flux=   55.4,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h07m00.20', '-47.54.21.7'),
          rad=17656,
          flux=   27.7,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h07m02.15', '-48.01.55.8'),
          rad=17898,
          flux=  176.2,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h07m02.72', '-44.43.21.5'),
          rad=14662,
          flux=  108.5,
          bmaj=  46.2,
          bmin=  0.0,
          bpa= -15.5)

addsource(t, dir=dm.direction('J2000', '12h07m08.26', '-44.16.32.6'),
          rad=14855,
          flux=   35.0,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h07m08.99', '-43.02.06.3'),
          rad=16371,
          flux=   17.1,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h07m11.80', '-42.36.15.5'),
          rad=17114,
          flux=   25.0,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h07m13.64', '-44.10.46.8'),
          rad=14875,
          flux=   14.3,
          bmaj=  29.0,
          bmin=  0.0,
          bpa=  58.2)

addsource(t, dir=dm.direction('J2000', '12h07m19.31', '-48.05.59.3'),
          rad=17901,
          flux=   82.0,
          bmaj=  17.0,
          bmin= 0.0,
          bpa= -42.8)

addsource(t, dir=dm.direction('J2000', '12h07m22.85', '-45.44.17.6'),
          rad=14531,
          flux=   18.3,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h07m23.77', '-42.23.13.1'),
          rad=17433,
          flux=   14.4,
          bmaj=  24.1,
          bmin=  0.0,
          bpa=  -2.1)

addsource(t, dir=dm.direction('J2000', '12h07m26.95', '-44.20.48.8'),
          rad=14608,
          flux=   22.6,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h07m29.24', '-45.32.30.5'),
          rad=14377,
          flux=   24.4,
          bmaj=  23.4,
          bmin=  0.0,
          bpa=  68.0)

addsource(t, dir=dm.direction('J2000', '12h07m32.35', '-46.19.14.4'),
          rad=14890,
          flux=   16.6,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h07m33.73', '-45.32.51.4'),
          rad=14332,
          flux=   24.4,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h07m47.21', '-45.14.00.0'),
          rad=14119,
          flux=  246.4,
          bmaj=  54.1,
          bmin=  0.0,
          bpa= -75.5)

addsource(t, dir=dm.direction('J2000', '12h07m52.60', '-45.43.16.2'),
          rad=14215,
          flux=   14.6,
          bmaj=  24.0,
          bmin= 0.0,
          bpa= -35.8)

addsource(t, dir=dm.direction('J2000', '12h07m57.65', '-47.29.01.9'),
          rad=16352,
          flux=   15.5,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h07m58.61', '-43.02.22.7'),
          rad=15884,
          flux=   35.6,
          bmaj=  17.5,
          bmin=  0.0,
          bpa= -83.5)

addsource(t, dir=dm.direction('J2000', '12h08m00.64', '-42.57.42.0'),
          rad=15999,
          flux=   16.5,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h08m00.74', '-43.39.20.6'),
          rad=14944,
          flux=   53.7,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h08m06.97', '-42.04.15.3'),
          rad=17721,
          flux=   43.5,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h08m07.49', '-47.31.14.2'),
          rad=16336,
          flux=   15.3,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h08m11.67', '-47.32.53.8'),
          rad=16352,
          flux=   17.9,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h08m12.26', '-44.33.31.5'),
          rad=14001,
          flux=   23.7,
          bmaj=  20.2,
          bmin=  0.0,
          bpa= -56.6)

addsource(t, dir=dm.direction('J2000', '12h08m13.39', '-44.03.43.3'),
          rad=14359,
          flux=   82.3,
          bmaj=  59.7,
          bmin= 0.0,
          bpa= -18.7)

addsource(t, dir=dm.direction('J2000', '12h08m13.56', '-42.43.04.5'),
          rad=16323,
          flux=   21.5,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h08m16.84', '-41.56.21.7'),
          rad=17933,
          flux=   11.3,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h08m20.70', '-43.29.01.5'),
          rad=14974,
          flux=   25.0,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h08m26.85', '-42.59.44.7'),
          rad=15689,
          flux=   26.4,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h08m38.12', '-47.11.23.1'),
          rad=15470,
          flux=   20.0,
          bmaj=  30.2,
          bmin=  0.0,
          bpa=  35.6)

addsource(t, dir=dm.direction('J2000', '12h08m41.00', '-45.49.42.7'),
          rad=13781,
          flux=   28.3,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h08m42.93', '-43.01.52.9'),
          rad=15473,
          flux=   11.7,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h08m47.04', '-43.37.52.7'),
          rad=14508,
          flux=   41.1,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h08m55.57', '-46.41.13.7'),
          rad=14526,
          flux=   64.5,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h08m55.89', '-46.39.46.9'),
          rad=14489,
          flux=   56.9,
          bmaj=  86.0,
          bmin=  65.0,
          bpa= -14.7)

addsource(t, dir=dm.direction('J2000', '12h08m56.37', '-46.43.09.0'),
          rad=14563,
          flux=   50.9,
          bmaj= 126.8,
          bmin=  55.6,
          bpa=  -9.9)

addsource(t, dir=dm.direction('J2000', '12h08m57.66', '-44.55.25.4'),
          rad=13389,
          flux=   18.6,
          bmaj=  51.0,
          bmin=  24.8,
          bpa=  35.5)

addsource(t, dir=dm.direction('J2000', '12h08m58.29', '-42.18.29.9'),
          rad=16751,
          flux=   24.1,
          bmaj=  20.1,
          bmin=  0.0,
          bpa=  33.7)

addsource(t, dir=dm.direction('J2000', '12h09m01.53', '-46.24.56.1'),
          rad=14118,
          flux=  110.4,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h09m09.25', '-46.10.18.6'),
          rad=13778,
          flux=   29.8,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h09m10.16', '-46.52.15.2'),
          rad=14659,
          flux=   25.9,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h09m13.52', '-43.39.18.7'),
          rad=14210,
          flux= 1979.9,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h09m14.90', '-43.59.58.0'),
          rad=13786,
          flux=   11.7,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h09m15.31', '-44.50.27.2'),
          rad=13221,
          flux=   16.3,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h09m15.61', '-46.18.46.5'),
          rad=13863,
          flux=   49.0,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h09m19.48', '-43.51.56.8'),
          rad=13887,
          flux=   11.5,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h09m19.74', '-46.31.06.2'),
          rad=14069,
          flux=   11.4,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h09m22.48', '-41.57.44.1'),
          rad=17323,
          flux=   11.6,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h09m27.20', '-48.31.24.8'),
          rad=17897,
          flux=   15.4,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h09m27.91', '-43.16.53.5'),
          rad=14619,
          flux=   13.2,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h09m28.76', '-44.26.37.4'),
          rad=13263,
          flux=   34.1,
          bmaj=  38.2,
          bmin=  0.0,
          bpa= -79.7)

addsource(t, dir=dm.direction('J2000', '12h09m30.24', '-43.41.04.2'),
          rad=14003,
          flux=   19.1,
          bmaj=  43.3,
          bmin=  25.1,
          bpa= -38.7)

addsource(t, dir=dm.direction('J2000', '12h09m30.91', '-48.25.02.0'),
          rad=17611,
          flux=   33.4,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h09m30.98', '-43.27.37.3'),
          rad=14311,
          flux=   37.4,
          bmaj=  17.4,
          bmin=  0.0,
          bpa=  23.6)

addsource(t, dir=dm.direction('J2000', '12h09m37.97', '-42.00.33.5'),
          rad=17083,
          flux=   11.9,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h09m39.27', '-45.00.32.8'),
          rad=12936,
          flux=   20.4,
          bmaj=  28.0,
          bmin=  26.6,
          bpa=  83.0)

addsource(t, dir=dm.direction('J2000', '12h09m42.77', '-48.34.29.0'),
          rad=17911,
          flux=   79.3,
          bmaj=  52.6,
          bmin=  0.0,
          bpa=  -8.9)

addsource(t, dir=dm.direction('J2000', '12h09m44.91', '-46.55.16.2'),
          rad=14419,
          flux=  144.9,
          bmaj=  24.7,
          bmin= 0.0,
          bpa= -30.9)

addsource(t, dir=dm.direction('J2000', '12h09m46.93', '-46.37.12.1'),
          rad=13946,
          flux=  118.2,
          bmaj=  36.6,
          bmin=  0.0,
          bpa= -18.9)

addsource(t, dir=dm.direction('J2000', '12h09m52.19', '-47.44.57.3'),
          rad=15922,
          flux=   21.9,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h09m57.57', '-46.18.06.1'),
          rad=13438,
          flux=  440.8,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h09m58.93', '-42.51.01.7'),
          rad=15091,
          flux=   27.9,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h10m00.42', '-44.24.40.7'),
          rad=12952,
          flux=  449.2,
          bmaj=  42.5,
          bmin=  0.0,
          bpa= -73.8)

addsource(t, dir=dm.direction('J2000', '12h10m12.32', '-44.03.31.6'),
          rad=13133,
          flux=   42.1,
          bmaj=  17.4,
          bmin= 0.0,
          bpa= -44.4)

addsource(t, dir=dm.direction('J2000', '12h10m15.56', '-45.31.31.1'),
          rad=12637,
          flux=   21.6,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h10m17.69', '-44.54.32.8'),
          rad=12545,
          flux=   15.1,
          bmaj=  22.5,
          bmin=  0.0,
          bpa= -21.7)

addsource(t, dir=dm.direction('J2000', '12h10m27.06', '-42.50.33.1'),
          rad=14847,
          flux=   15.1,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h10m27.93', '-47.50.03.1'),
          rad=15822,
          flux=   26.1,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h10m28.10', '-48.12.56.9'),
          rad=16704,
          flux=   37.6,
          bmaj=  15.8,
          bmin=  0.0,
          bpa= -32.9)

addsource(t, dir=dm.direction('J2000', '12h10m32.98', '-44.47.30.1'),
          rad=12414,
          flux=   34.4,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h10m34.77', '-47.40.30.8'),
          rad=15420,
          flux=   35.7,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h10m36.72', '-46.06.00.7'),
          rad=12835,
          flux=   40.5,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h10m40.17', '-45.28.53.5'),
          rad=12363,
          flux=   46.2,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h10m43.55', '-48.34.24.3'),
          rad=17481,
          flux=   96.6,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h10m47.67', '-41.56.00.1'),
          rad=16687,
          flux=   46.1,
          bmaj=  24.6,
          bmin=  0.0,
          bpa=  31.3)

addsource(t, dir=dm.direction('J2000', '12h10m52.08', '-44.41.40.2'),
          rad=12248,
          flux=   30.2,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h10m55.38', '-42.33.03.4'),
          rad=15195,
          flux=   70.8,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h10m57.97', '-42.29.31.2'),
          rad=15300,
          flux=   21.6,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h11m02.17', '-42.17.07.2'),
          rad=15730,
          flux=   14.5,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h11m02.79', '-45.39.12.6'),
          rad=12212,
          flux=  231.9,
          bmaj=  39.3,
          bmin=  0.0,
          bpa= -43.0)

addsource(t, dir=dm.direction('J2000', '12h11m03.52', '-44.00.18.9'),
          rad=12665,
          flux=  128.6,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h11m05.66', '-44.35.53.0'),
          rad=12151,
          flux=   32.8,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h11m07.12', '-47.09.02.5'),
          rad=14089,
          flux=   22.3,
          bmaj=  18.1,
          bmin=  0.0,
          bpa= -83.0)

addsource(t, dir=dm.direction('J2000', '12h11m11.75', '-45.53.03.1'),
          rad=12284,
          flux=   11.5,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h11m16.12', '-45.45.47.0'),
          rad=12146,
          flux=   30.2,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h11m18.51', '-42.34.26.0'),
          rad=14942,
          flux=  334.6,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h11m21.92', '-42.18.00.5'),
          rad=15528,
          flux=   16.3,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h11m27.00', '-47.03.32.8'),
          rad=13745,
          flux=   60.0,
          bmaj=  15.9,
          bmin=  0.0,
          bpa= -38.3)

addsource(t, dir=dm.direction('J2000', '12h11m28.28', '-48.31.52.4'),
          rad=17064,
          flux=   24.4,
          bmaj=  18.4,
          bmin= 0.0,
          bpa= -36.0)

addsource(t, dir=dm.direction('J2000', '12h11m37.04', '-45.56.09.6'),
          rad=12073,
          flux=   18.7,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h11m48.07', '-45.58.36.0'),
          rad=12000,
          flux=   17.6,
          bmaj=  32.3,
          bmin=  0.0,
          bpa=  76.3)

addsource(t, dir=dm.direction('J2000', '12h11m50.00', '-47.12.36.1'),
          rad=13835,
          flux=   67.0,
          bmaj=  30.3,
          bmin= 0.0,
          bpa=  51.8)

addsource(t, dir=dm.direction('J2000', '12h11m50.04', '-43.52.03.2'),
          rad=12356,
          flux=   10.8,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h11m51.06', '-48.49.15.2'),
          rad=17684,
          flux=  135.2,
          bmaj=  30.9,
          bmin=  0.0,
          bpa=  39.7)

addsource(t, dir=dm.direction('J2000', '12h11m51.98', '-42.54.39.9'),
          rad=13935,
          flux=   27.8,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h11m52.57', '-45.44.14.9'),
          rad=11755,
          flux=   21.0,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h11m53.84', '-41.49.24.4'),
          rad=16435,
          flux=   81.1,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h11m58.23', '-47.08.59.2'),
          rad=13647,
          flux=   13.3,
          bmaj=  25.3,
          bmin=  0.0,
          bpa= -46.4)

addsource(t, dir=dm.direction('J2000', '12h11m59.38', '-47.28.42.1'),
          rad=14312,
          flux=   12.5,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h12m00.07', '-47.26.30.0'),
          rad=14228,
          flux=   27.5,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h12m02.66', '-46.56.14.0'),
          rad=13209,
          flux=  100.2,
          bmaj=  19.5,
          bmin=  0.0,
          bpa=  48.3)

addsource(t, dir=dm.direction('J2000', '12h12m03.33', '-46.08.12.0'),
          rad=12015,
          flux=   11.1,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h12m03.89', '-41.56.03.1'),
          rad=16074,
          flux=   57.4,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h12m04.86', '-47.00.00.1'),
          rad=13305,
          flux=   21.2,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h12m06.66', '-45.53.18.7'),
          rad=11732,
          flux=   16.0,
          bmaj=  22.5,
          bmin=  0.0,
          bpa= -46.6)

addsource(t, dir=dm.direction('J2000', '12h12m10.18', '-42.23.00.1'),
          rad=14929,
          flux=   21.8,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h12m14.59', '-44.02.04.6'),
          rad=11905,
          flux=   39.1,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h12m18.08', '-42.15.47.2'),
          rad=15148,
          flux=   23.2,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h12m20.40', '-48.31.24.8'),
          rad=16693,
          flux=   57.5,
          bmaj=  34.8,
          bmin=  0.0,
          bpa= -88.7)

addsource(t, dir=dm.direction('J2000', '12h12m22.69', '-48.19.10.9'),
          rad=16143,
          flux=   12.5,
          bmaj=  25.2,
          bmin= 0.0,
          bpa= -29.7)

addsource(t, dir=dm.direction('J2000', '12h12m22.75', '-44.02.49.9'),
          rad=11807,
          flux=   36.8,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h12m32.25', '-42.58.22.0'),
          rad=13446,
          flux=   46.8,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h12m37.79', '-47.50.04.5'),
          rad=14825,
          flux=   58.2,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h12m39.73', '-42.23.50.5'),
          rad=14649,
          flux=   24.0,
          bmaj=  23.5,
          bmin=  0.0,
          bpa= -22.1)

addsource(t, dir=dm.direction('J2000', '12h12m39.77', '-44.21.52.5'),
          rad=11320,
          flux=   34.2,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h12m48.39', '-42.20.44.9'),
          rad=14700,
          flux=  191.5,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h12m48.78', '-45.43.40.7'),
          rad=11173,
          flux=   16.8,
          bmaj=  21.2,
          bmin=  0.0,
          bpa=  24.7)

addsource(t, dir=dm.direction('J2000', '12h12m52.95', '-44.04.56.4'),
          rad=11459,
          flux=   14.4,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h12m57.77', '-48.48.24.8'),
          rad=17224,
          flux=  204.4,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h13m02.28', '-41.10.14.8'),
          rad=17702,
          flux=   29.4,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h13m03.24', '-43.19.41.2'),
          rad=12477,
          flux=  142.2,
          bmaj=  15.2,
          bmin=  0.0,
          bpa= -85.2)

addsource(t, dir=dm.direction('J2000', '12h13m07.77', '-42.04.07.3'),
          rad=15231,
          flux=   37.4,
          bmaj=  22.3,
          bmin= 0.0,
          bpa=  11.0)

addsource(t, dir=dm.direction('J2000', '12h13m08.98', '-42.02.33.3'),
          rad=15288,
          flux=   42.8,
          bmaj=  19.0,
          bmin=  0.0,
          bpa= -39.2)

addsource(t, dir=dm.direction('J2000', '12h13m11.37', '-42.33.10.9'),
          rad=14019,
          flux=   35.4,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h13m15.25', '-45.05.39.1'),
          rad=10648,
          flux=   14.0,
          bmaj=  24.9,
          bmin=  0.0,
          bpa= -40.6)

addsource(t, dir=dm.direction('J2000', '12h13m26.05', '-43.55.34.8'),
          rad=11314,
          flux=   10.7,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h13m26.07', '-45.17.03.1'),
          rad=10560,
          flux=  120.5,
          bmaj=  14.6,
          bmin=  0.0,
          bpa= -46.8)

addsource(t, dir=dm.direction('J2000', '12h13m27.44', '-48.23.13.3'),
          rad=15884,
          flux=  256.2,
          bmaj=  18.9,
          bmin= 0.0,
          bpa=  88.7)

addsource(t, dir=dm.direction('J2000', '12h13m30.04', '-45.58.59.3'),
          rad=10988,
          flux=   15.2,
          bmaj=  35.5,
          bmin=  0.0,
          bpa=   2.8)

addsource(t, dir=dm.direction('J2000', '12h13m30.69', '-43.25.12.9'),
          rad=12053,
          flux=   73.3,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h13m32.04', '-48.23.14.5'),
          rad=15855,
          flux=  167.8,
          bmaj=  19.4,
          bmin=  0.0,
          bpa=  84.9)

addsource(t, dir=dm.direction('J2000', '12h13m33.91', '-41.37.09.4'),
          rad=16226,
          flux=   17.1,
          bmaj=  40.7,
          bmin= 0.0,
          bpa= -86.9)

addsource(t, dir=dm.direction('J2000', '12h13m40.84', '-43.42.31.3'),
          rad=11477,
          flux=   26.4,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h13m40.86', '-43.30.43.0'),
          rad=11796,
          flux=   12.9,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h13m41.05', '-48.23.25.2'),
          rad=15804,
          flux=   12.3,
          bmaj=  31.8,
          bmin=  0.0,
          bpa=   1.3)

addsource(t, dir=dm.direction('J2000', '12h13m45.59', '-43.34.58.8'),
          rad=11632,
          flux=   85.1,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h13m47.71', '-45.30.39.8'),
          rad=10424,
          flux=   11.9,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h13m48.00', '-46.28.46.5'),
          rad=11477,
          flux=   31.6,
          bmaj=  35.6,
          bmin=  0.0,
          bpa=  46.3)

addsource(t, dir=dm.direction('J2000', '12h13m56.15', '-43.03.20.8'),
          rad=12522,
          flux=   37.6,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h13m58.54', '-48.43.30.2'),
          rad=16622,
          flux=   50.0,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h13m58.89', '-41.17.09.6'),
          rad=16988,
          flux=   17.1,
          bmaj=  21.1,
          bmin=  0.0,
          bpa= -62.5)

addsource(t, dir=dm.direction('J2000', '12h14m00.62', '-45.23.49.5'),
          rad=10235,
          flux=   30.0,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h14m08.95', '-47.03.49.1'),
          rad=12371,
          flux=   14.8,
          bmaj=  23.8,
          bmin=  0.0,
          bpa=  23.2)

addsource(t, dir=dm.direction('J2000', '12h14m11.09', '-48.49.31.0'),
          rad=16834,
          flux=   19.9,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h14m14.58', '-41.10.26.9'),
          rad=17208,
          flux=   28.0,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h14m15.22', '-45.40.04.7'),
          rad=10243,
          flux=   92.7,
          bmaj=  62.7,
          bmin=  0.0,
          bpa=  69.7)

addsource(t, dir=dm.direction('J2000', '12h14m17.93', '-49.08.51.9'),
          rad=17737,
          flux=   23.9,
          bmaj=  19.4,
          bmin= 0.0,
          bpa=  39.6)

addsource(t, dir=dm.direction('J2000', '12h14m21.36', '-43.54.42.8'),
          rad=10780,
          flux=   57.9,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h14m22.04', '-48.46.53.2'),
          rad=16643,
          flux=   10.9,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h14m26.36', '-48.19.55.2'),
          rad=15352,
          flux=  105.3,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h14m32.97', '-41.21.32.4'),
          rad=16550,
          flux=   57.5,
          bmaj=  27.3,
          bmin=  0.0,
          bpa= -70.8)

addsource(t, dir=dm.direction('J2000', '12h14m33.28', '-41.12.36.4'),
          rad=16981,
          flux=   69.8,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h14m35.00', '-42.32.10.0'),
          rad=13370,
          flux=   57.8,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h14m36.90', '-41.11.40.1'),
          rad=17003,
          flux=   44.5,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h14m37.07', '-48.08.59.3'),
          rad=14787,
          flux=   13.8,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h14m43.16', '-43.31.03.1'),
          rad=11195,
          flux=   27.6,
          bmaj=  31.0,
          bmin=  0.0,
          bpa= -59.3)

addsource(t, dir=dm.direction('J2000', '12h14m49.86', '-42.19.36.6'),
          rad=13776,
          flux=  113.8,
          bmaj=  38.1,
          bmin=  18.7,
          bpa= -18.6)

addsource(t, dir=dm.direction('J2000', '12h14m50.40', '-44.53.17.6'),
          rad= 9661,
          flux=   12.2,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h14m52.16', '-42.18.10.7'),
          rad=13820,
          flux=  122.9,
          bmaj=  46.1,
          bmin=  18.0,
          bpa= -23.7)

addsource(t, dir=dm.direction('J2000', '12h14m52.30', '-45.27.47.6'),
          rad= 9728,
          flux=   48.4,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h15m00.45', '-48.21.39.3'),
          rad=15217,
          flux=   11.8,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h15m02.26', '-47.46.40.6'),
          rad=13635,
          flux=   20.9,
          bmaj=  46.5,
          bmin=  0.0,
          bpa=  47.5)

addsource(t, dir=dm.direction('J2000', '12h15m09.99', '-46.27.15.9'),
          rad=10683,
          flux=  271.5,
          bmaj=  16.6,
          bmin= 0.0,
          bpa=  25.5)

addsource(t, dir=dm.direction('J2000', '12h15m12.79', '-45.53.29.4'),
          rad= 9868,
          flux=   26.0,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h15m16.03', '-45.21.27.1'),
          rad= 9431,
          flux=   21.7,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h15m17.50', '-41.35.58.1'),
          rad=15560,
          flux=   37.5,
          bmaj=  40.1,
          bmin=  0.0,
          bpa= -22.4)

addsource(t, dir=dm.direction('J2000', '12h15m19.31', '-47.21.43.8'),
          rad=12479,
          flux=   32.4,
          bmaj=  20.0,
          bmin= 0.0,
          bpa=  58.4)

addsource(t, dir=dm.direction('J2000', '12h15m21.04', '-45.33.34.1'),
          rad= 9489,
          flux=   50.6,
          bmaj=  23.1,
          bmin=  0.0,
          bpa=  33.4)

addsource(t, dir=dm.direction('J2000', '12h15m21.72', '-47.17.30.4'),
          rad=12294,
          flux=   17.4,
          bmaj=  30.1,
          bmin=  21.8,
          bpa= -37.7)

addsource(t, dir=dm.direction('J2000', '12h15m22.65', '-41.46.18.7'),
          rad=15034,
          flux=   15.6,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h15m23.51', '-47.02.00.8'),
          rad=11694,
          flux=   14.0,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h15m24.68', '-44.34.11.3'),
          rad= 9443,
          flux=   33.4,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h15m27.20', '-41.34.33.6'),
          rad=15562,
          flux=   12.9,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h15m30.96', '-45.09.48.7'),
          rad= 9219,
          flux=   24.0,
          bmaj=  23.9,
          bmin=  0.0,
          bpa= -42.6)

addsource(t, dir=dm.direction('J2000', '12h15m35.98', '-47.46.11.6'),
          rad=13380,
          flux=   11.5,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h15m36.88', '-42.11.47.1'),
          rad=13763,
          flux=   38.5,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h15m38.76', '-46.53.08.8'),
          rad=11252,
          flux=   13.3,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h15m40.95', '-45.17.12.5'),
          rad= 9143,
          flux=   50.4,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h15m41.22', '-48.51.07.2'),
          rad=16401,
          flux=   12.9,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h15m44.34', '-45.29.12.8'),
          rad= 9202,
          flux=   65.8,
          bmaj=  16.6,
          bmin=  0.0,
          bpa= -12.3)

addsource(t, dir=dm.direction('J2000', '12h15m46.96', '-40.50.29.9'),
          rad=17635,
          flux=   12.6,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h15m47.44', '-41.29.06.8'),
          rad=15694,
          flux=   31.0,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h15m49.05', '-46.03.52.5'),
          rad= 9723,
          flux=   31.3,
          bmaj=  17.4,
          bmin=  0.0,
          bpa= -61.8)

addsource(t, dir=dm.direction('J2000', '12h15m50.65', '-47.05.26.1'),
          rad=11604,
          flux=  120.3,
          bmaj= 111.9,
          bmin=  30.2,
          bpa=  50.3)

addsource(t, dir=dm.direction('J2000', '12h15m50.95', '-48.00.32.5'),
          rad=13921,
          flux=   15.2,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h15m56.88', '-48.52.34.3'),
          rad=16388,
          flux=   57.3,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h16m00.67', '-42.01.27.7'),
          rad=14062,
          flux=   57.5,
          bmaj=  31.6,
          bmin=  0.0,
          bpa= -38.7)

addsource(t, dir=dm.direction('J2000', '12h16m04.57', '-43.23.35.7'),
          rad=10679,
          flux=   12.4,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h16m04.72', '-43.40.20.9'),
          rad=10151,
          flux=   15.5,
          bmaj=  25.8,
          bmin=  0.0,
          bpa= -89.6)

addsource(t, dir=dm.direction('J2000', '12h16m04.83', '-42.00.41.0'),
          rad=14069,
          flux=  177.3,
          bmaj= 127.2,
          bmin=  41.7,
          bpa= -36.4)

addsource(t, dir=dm.direction('J2000', '12h16m13.55', '-41.39.20.5'),
          rad=15026,
          flux=   14.7,
          bmaj=  26.8,
          bmin= 0.0,
          bpa=  85.1)

addsource(t, dir=dm.direction('J2000', '12h16m15.52', '-42.37.44.0'),
          rad=12340,
          flux=  212.4,
          bmaj=  25.7,
          bmin=  0.0,
          bpa=  -0.5)

addsource(t, dir=dm.direction('J2000', '12h16m16.07', '-41.58.02.3'),
          rad=14114,
          flux=  198.6,
          bmaj=  66.9,
          bmin=  40.5,
          bpa= -41.2)

addsource(t, dir=dm.direction('J2000', '12h16m17.20', '-47.13.37.5'),
          rad=11717,
          flux=   19.4,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h16m21.62', '-41.50.51.6'),
          rad=14417,
          flux=   21.7,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h16m23.29', '-44.09.01.8'),
          rad= 9243,
          flux=   10.8,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h16m23.82', '-44.48.41.5'),
          rad= 8694,
          flux=   21.4,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h16m23.98', '-49.21.24.3'),
          rad=17728,
          flux=   11.2,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h16m23.99', '-43.01.47.6'),
          rad=11297,
          flux=   13.7,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h16m26.33', '-47.01.32.4'),
          rad=11174,
          flux=   22.2,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h16m26.58', '-42.23.32.9'),
          rad=12871,
          flux=   41.6,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h16m40.80', '-46.30.13.5'),
          rad= 9958,
          flux=   27.1,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h16m41.81', '-41.10.54.5'),
          rad=16270,
          flux=   16.3,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h16m43.19', '-46.07.44.6'),
          rad= 9298,
          flux=   45.5,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h16m44.95', '-45.14.12.7'),
          rad= 8455,
          flux=  131.3,
          bmaj=  28.9,
          bmin=  0.0,
          bpa= -83.2)

addsource(t, dir=dm.direction('J2000', '12h16m45.70', '-41.39.53.9'),
          rad=14792,
          flux=  152.9,
          bmaj=  34.3,
          bmin= 0.0,
          bpa=  31.2)

addsource(t, dir=dm.direction('J2000', '12h16m46.75', '-44.46.53.0'),
          rad= 8463,
          flux=   12.2,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h16m48.74', '-43.31.06.3'),
          rad=10030,
          flux=   12.1,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h16m51.96', '-41.38.42.9'),
          rad=14810,
          flux=  129.2,
          bmaj=  20.8,
          bmin=  0.0,
          bpa= -42.4)

addsource(t, dir=dm.direction('J2000', '12h16m52.86', '-47.43.50.6'),
          rad=12757,
          flux=   19.0,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h16m53.32', '-40.41.52.0'),
          rad=17715,
          flux=   34.6,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h16m53.86', '-45.44.35.3'),
          rad= 8702,
          flux=   12.2,
          bmaj=  28.4,
          bmin=  0.0,
          bpa=  77.3)

addsource(t, dir=dm.direction('J2000', '12h17m00.21', '-47.15.13.1'),
          rad=11461,
          flux=   90.6,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h17m01.23', '-46.55.50.2'),
          rad=10682,
          flux=   16.2,
          bmaj=  26.9,
          bmin=  0.0,
          bpa= -39.6)

addsource(t, dir=dm.direction('J2000', '12h17m11.68', '-43.37.18.3'),
          rad= 9620,
          flux=   55.4,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h17m14.82', '-48.58.19.6'),
          rad=16282,
          flux=   44.9,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h17m15.30', '-47.12.14.5'),
          rad=11227,
          flux=   94.3,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h17m18.58', '-42.37.42.9'),
          rad=11858,
          flux=   13.5,
          bmaj=  28.9,
          bmin=  0.0,
          bpa= -24.9)

addsource(t, dir=dm.direction('J2000', '12h17m20.60', '-42.12.11.8'),
          rad=13005,
          flux=   58.6,
          bmaj=  30.6,
          bmin= 0.0,
          bpa= -70.8)

addsource(t, dir=dm.direction('J2000', '12h17m21.06', '-45.19.48.3'),
          rad= 8111,
          flux=  251.0,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h17m21.32', '-41.08.19.5'),
          rad=16177,
          flux=   11.7,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h17m23.53', '-45.54.08.2'),
          rad= 8594,
          flux=   13.0,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h17m25.94', '-49.06.02.8'),
          rad=16632,
          flux=   25.9,
          bmaj=  24.3,
          bmin=  0.0,
          bpa= -46.1)

addsource(t, dir=dm.direction('J2000', '12h17m26.57', '-40.34.34.1'),
          rad=17930,
          flux=   54.9,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h17m27.74', '-42.22.24.4'),
          rad=12480,
          flux=   45.7,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h17m28.27', '-45.37.27.3'),
          rad= 8239,
          flux=   15.1,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h17m31.82', '-48.59.54.2'),
          rad=16282,
          flux=   44.1,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h17m36.94', '-43.09.30.3'),
          rad=10390,
          flux=   23.1,
          bmaj=  25.0,
          bmin=  0.0,
          bpa=  73.0)

addsource(t, dir=dm.direction('J2000', '12h17m42.71', '-46.50.39.5'),
          rad=10157,
          flux=   13.9,
          bmaj=  38.4,
          bmin= 0.0,
          bpa= -81.1)

addsource(t, dir=dm.direction('J2000', '12h17m43.64', '-49.22.35.1'),
          rad=17427,
          flux=   15.9,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h17m44.67', '-41.04.48.4'),
          rad=16232,
          flux=   47.8,
          bmaj=  52.4,
          bmin=  16.9,
          bpa=  84.7)

addsource(t, dir=dm.direction('J2000', '12h17m48.48', '-46.34.53.2'),
          rad= 9532,
          flux=   15.5,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h17m49.55', '-49.02.06.9'),
          rad=16315,
          flux=   91.3,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h17m53.73', '-45.42.31.9'),
          rad= 8067,
          flux=  133.7,
          bmaj=  45.8,
          bmin=  30.6,
          bpa= -89.4)

addsource(t, dir=dm.direction('J2000', '12h17m54.15', '-49.08.27.4'),
          rad=16629,
          flux=   34.8,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h17m59.25', '-42.32.48.0'),
          rad=11779,
          flux=   15.7,
          bmaj=  28.0,
          bmin=  0.0,
          bpa=  76.3)

addsource(t, dir=dm.direction('J2000', '12h18m02.59', '-40.51.53.6'),
          rad=16821,
          flux=   14.5,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h18m05.76', '-42.36.36.8'),
          rad=11559,
          flux=   22.5,
          bmaj=  19.6,
          bmin=  0.0,
          bpa= -74.1)

addsource(t, dir=dm.direction('J2000', '12h18m05.89', '-46.41.19.8'),
          rad= 9620,
          flux=   28.9,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h18m13.93', '-43.27.23.6'),
          rad= 9401,
          flux=  121.3,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h18m13.97', '-49.07.37.6'),
          rad=16495,
          flux=   11.6,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h18m15.53', '-42.54.04.2'),
          rad=10716,
          flux=   18.2,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h18m20.66', '-45.43.34.3'),
          rad= 7818,
          flux=   33.6,
          bmaj=  21.9,
          bmin=  0.0,
          bpa=  51.1)

addsource(t, dir=dm.direction('J2000', '12h18m25.15', '-41.09.54.9'),
          rad=15747,
          flux=   14.4,
          bmaj=  25.3,
          bmin= 0.0,
          bpa=  19.1)

addsource(t, dir=dm.direction('J2000', '12h18m26.69', '-45.45.04.4'),
          rad= 7787,
          flux=   40.6,
          bmaj=  24.2,
          bmin=  0.0,
          bpa=  23.4)

addsource(t, dir=dm.direction('J2000', '12h18m28.11', '-43.04.25.1'),
          rad=10181,
          flux=   12.8,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h18m29.61', '-41.56.45.3'),
          rad=13306,
          flux=   18.7,
          bmaj=  25.8,
          bmin=  0.0,
          bpa= -45.0)

addsource(t, dir=dm.direction('J2000', '12h18m35.83', '-43.54.26.3'),
          rad= 8312,
          flux=   23.8,
          bmaj=  22.0,
          bmin= 0.0,
          bpa=  20.9)

addsource(t, dir=dm.direction('J2000', '12h18m45.64', '-47.19.03.1'),
          rad=10888,
          flux=  202.4,
          bmaj= 113.1,
          bmin=  82.4,
          bpa=  45.7)

addsource(t, dir=dm.direction('J2000', '12h18m46.40', '-47.14.26.6'),
          rad=10676,
          flux=   16.9,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h18m48.13', '-41.51.29.3'),
          rad=13461,
          flux=  100.3,
          bmaj=  36.1,
          bmin=  0.0,
          bpa=  48.9)

addsource(t, dir=dm.direction('J2000', '12h18m48.44', '-44.02.00.3'),
          rad= 7979,
          flux=   22.8,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h18m48.50', '-49.35.27.5'),
          rad=17858,
          flux=   23.9,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h18m55.10', '-47.18.47.7'),
          rad=10814,
          flux=  238.9,
          bmaj= 123.2,
          bmin=  51.5,
          bpa= -46.2)

addsource(t, dir=dm.direction('J2000', '12h18m55.39', '-47.43.18.9'),
          rad=11964,
          flux=   31.4,
          bmaj=  0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h18m57.38', '-45.37.00.4'),
          rad= 7332,
          flux=   28.7,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h18m59.30', '-48.09.57.1'),
          rad=13266,
          flux=   42.3,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h19m00.87', '-46.25.09.6'),
          rad= 8586,
          flux=   38.2,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h19m04.61', '-43.38.08.5'),
          rad= 8576,
          flux=   11.3,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h19m04.87', '-46.30.19.1'),
          rad= 8736,
          flux=   15.9,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h19m06.50', '-48.29.57.0'),
          rad=14262,
          flux=  618.5,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h19m07.33', '-47.17.10.7'),
          rad=10659,
          flux=   43.5,
          bmaj=  96.8,
          bmin=  18.2,
          bpa=  74.0)

addsource(t, dir=dm.direction('J2000', '12h19m07.82', '-42.54.35.3'),
          rad=10301,
          flux=   12.6,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h19m17.10', '-46.06.34.7'),
          rad= 7843,
          flux=   25.7,
          bmaj=  18.0,
          bmin= 0.0,
          bpa=  59.7)

addsource(t, dir=dm.direction('J2000', '12h19m18.29', '-44.54.15.9'),
          rad= 6819,
          flux=   36.7,
          bmaj=  18.4,
          bmin= 0.0,
          bpa= -49.1)

addsource(t, dir=dm.direction('J2000', '12h19m19.95', '-45.40.22.8'),
          rad= 7169,
          flux=   17.1,
          bmaj=  33.9,
          bmin= 0.0,
          bpa=  19.8)

addsource(t, dir=dm.direction('J2000', '12h19m20.10', '-48.00.01.7'),
          rad=12652,
          flux=   21.2,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h19m24.19', '-44.47.13.1'),
          rad= 6798,
          flux=   66.7,
          bmaj=  16.1,
          bmin= 0.0,
          bpa=  88.3)

addsource(t, dir=dm.direction('J2000', '12h19m24.21', '-45.38.16.5'),
          rad= 7086,
          flux=   19.3,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h19m31.83', '-46.44.50.0'),
          rad= 9085,
          flux=   20.3,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h19m32.27', '-47.18.07.5'),
          rad=10541,
          flux=   14.4,
          bmaj=  31.5,
          bmin= 0.0,
          bpa=   9.5)

addsource(t, dir=dm.direction('J2000', '12h19m39.89', '-45.36.04.0'),
          rad= 6890,
          flux=   70.2,
          bmaj=  15.2,
          bmin= 0.0,
          bpa=   7.8)

addsource(t, dir=dm.direction('J2000', '12h19m42.29', '-45.56.45.6'),
          rad= 7334,
          flux=   16.6,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h19m47.26', '-48.00.47.3'),
          rad=12548,
          flux=   11.3,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h19m48.26', '-42.01.09.0'),
          rad=12618,
          flux=  110.4,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h19m52.33', '-46.54.26.2'),
          rad= 9339,
          flux=   22.7,
          bmaj=  24.7,
          bmin= 0.0,
          bpa= -52.5)

addsource(t, dir=dm.direction('J2000', '12h19m53.83', '-42.29.51.2'),
          rad=11142,
          flux=  302.8,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h19m57.28', '-48.23.37.8'),
          rad=13689,
          flux=   52.1,
          bmaj=  77.4,
          bmin=  20.5,
          bpa=  11.1)

addsource(t, dir=dm.direction('J2000', '12h19m57.41', '-40.34.30.3'),
          rad=17233,
          flux=   63.2,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h19m58.40', '-46.25.32.7'),
          rad= 8124,
          flux=   27.4,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h19m59.80', '-45.47.45.0'),
          rad= 6939,
          flux=   13.1,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h19m59.96', '-48.52.26.9'),
          rad=15225,
          flux=   10.9,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h20m02.38', '-40.25.40.4'),
          rad=17704,
          flux=   31.7,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h20m03.21', '-40.26.35.4'),
          rad=17650,
          flux=   35.3,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h20m04.25', '-47.21.44.1'),
          rad=10511,
          flux=   10.6,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h20m05.98', '-42.16.41.7'),
          rad=11722,
          flux=   15.6,
          bmaj=  39.9,
          bmin=  27.6,
          bpa= -65.4)

addsource(t, dir=dm.direction('J2000', '12h20m06.36', '-40.42.30.0'),
          rad=16752,
          flux=   72.3,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h20m10.87', '-48.05.26.1'),
          rad=12668,
          flux=   47.0,
          bmaj=  70.6,
          bmin= 0.0,
          bpa=   1.2)

addsource(t, dir=dm.direction('J2000', '12h20m13.29', '-47.59.53.9'),
          rad=12369,
          flux=   51.8,
          bmaj=  30.5,
          bmin= 0.0,
          bpa=  50.2)

addsource(t, dir=dm.direction('J2000', '12h20m13.74', '-41.40.14.3'),
          rad=13574,
          flux=   68.0,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h20m15.84', '-47.39.25.4'),
          rad=11311,
          flux=   16.1,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h20m17.10', '-46.03.32.7'),
          rad= 7213,
          flux=   33.6,
          bmaj=  70.8,
          bmin= 0.0,
          bpa=   5.7)

addsource(t, dir=dm.direction('J2000', '12h20m18.13', '-48.46.42.2'),
          rad=14838,
          flux=   35.9,
          bmaj=  17.8,
          bmin= 0.0,
          bpa=  71.9)

addsource(t, dir=dm.direction('J2000', '12h20m19.23', '-43.18.32.0'),
          rad= 8722,
          flux=   40.8,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h20m23.79', '-41.58.33.8'),
          rad=12554,
          flux=   17.2,
          bmaj=  34.8,
          bmin= 0.0,
          bpa= -13.4)

addsource(t, dir=dm.direction('J2000', '12h20m24.55', '-43.10.28.1'),
          rad= 9031,
          flux=   28.7,
          bmaj=  42.4,
          bmin= 0.0,
          bpa= -52.2)

addsource(t, dir=dm.direction('J2000', '12h20m25.06', '-44.45.41.1'),
          rad= 6170,
          flux=   37.7,
          bmaj=  54.2,
          bmin= 0.0,
          bpa=  16.8)

addsource(t, dir=dm.direction('J2000', '12h20m27.61', '-47.53.24.8'),
          rad=11961,
          flux=   24.6,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h20m29.80', '-42.19.10.9'),
          rad=11455,
          flux=   17.7,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h20m34.91', '-41.32.20.2'),
          rad=13892,
          flux=   28.0,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h20m35.91', '-43.04.13.5'),
          rad= 9229,
          flux=   34.6,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h20m43.68', '-47.09.30.4'),
          rad= 9685,
          flux=   47.7,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h20m46.06', '-49.45.10.0'),
          rad=17986,
          flux=   14.9,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h20m46.52', '-41.07.46.4'),
          rad=15181,
          flux=  150.3,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h20m49.29', '-44.02.41.8'),
          rad= 6818,
          flux=   10.8,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h20m50.17', '-41.50.36.5'),
          rad=12836,
          flux=   30.1,
          bmaj=  42.4,
          bmin= 0.0,
          bpa= -46.4)

addsource(t, dir=dm.direction('J2000', '12h20m53.11', '-45.36.42.2'),
          rad= 6174,
          flux=   33.7,
          bmaj=  20.0,
          bmin= 0.0,
          bpa=  69.9)

addsource(t, dir=dm.direction('J2000', '12h20m56.16', '-48.42.34.9'),
          rad=14460,
          flux=   11.4,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h20m58.02', '-45.04.22.8'),
          rad= 5750,
          flux=   11.5,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h20m59.32', '-41.49.12.1'),
          rad=12865,
          flux=   46.5,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h21m02.28', '-45.06.00.7'),
          rad= 5709,
          flux=   17.5,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h21m03.07', '-40.26.01.3'),
          rad=17449,
          flux=   15.2,
          bmaj=  24.5,
          bmin= 0.0,
          bpa=  76.2)

addsource(t, dir=dm.direction('J2000', '12h21m03.39', '-41.40.24.2'),
          rad=13320,
          flux=   17.3,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h21m06.05', '-45.48.45.5'),
          rad= 6337,
          flux=  123.1,
          bmaj=  31.7,
          bmin= 0.0,
          bpa=  40.2)

addsource(t, dir=dm.direction('J2000', '12h21m07.79', '-47.42.29.7'),
          rad=11192,
          flux=   41.9,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h21m10.92', '-48.42.39.8'),
          rad=14407,
          flux=  275.3,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h21m12.63', '-41.52.18.9'),
          rad=12632,
          flux=  505.1,
          bmaj=  46.2,
          bmin=  17.8,
          bpa= -44.2)

addsource(t, dir=dm.direction('J2000', '12h21m12.99', '-44.10.10.3'),
          rad= 6373,
          flux=   55.2,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h21m13.70', '-41.14.32.2'),
          rad=14690,
          flux=  164.0,
          bmaj=  19.0,
          bmin= 0.0,
          bpa=  38.3)

addsource(t, dir=dm.direction('J2000', '12h21m20.34', '-47.50.14.3'),
          rad=11534,
          flux=   11.1,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h21m27.03', '-44.11.29.8'),
          rad= 6203,
          flux=  557.5,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h21m27.34', '-47.53.28.6'),
          rad=11672,
          flux=   35.2,
          bmaj=  50.4,
          bmin= 0.0,
          bpa=  43.7)

addsource(t, dir=dm.direction('J2000', '12h21m27.35', '-49.18.26.7'),
          rad=16346,
          flux=   28.4,
          bmaj=  61.2,
          bmin= 0.0,
          bpa= -61.6)

addsource(t, dir=dm.direction('J2000', '12h21m28.46', '-47.50.59.8'),
          rad=11536,
          flux=   16.2,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h21m30.34', '-43.23.37.6'),
          rad= 7965,
          flux=   16.0,
          bmaj=  22.6,
          bmin= 0.0,
          bpa= -45.5)

addsource(t, dir=dm.direction('J2000', '12h21m30.79', '-41.29.27.2'),
          rad=13792,
          flux=   64.9,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h21m35.62', '-47.32.46.7'),
          rad=10548,
          flux=   35.4,
          bmaj=  37.7,
          bmin= 0.0,
          bpa= -37.9)

addsource(t, dir=dm.direction('J2000', '12h21m37.30', '-49.41.23.6'),
          rad=17616,
          flux=   22.4,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h21m40.47', '-45.45.15.0'),
          rad= 5921,
          flux=   17.9,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h21m40.61', '-41.46.15.0'),
          rad=12827,
          flux=   28.3,
          bmaj=  23.2,
          bmin= 0.0,
          bpa=  34.1)

addsource(t, dir=dm.direction('J2000', '12h21m42.25', '-40.15.18.9'),
          rad=17918,
          flux=   31.4,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h21m42.50', '-43.53.00.2'),
          rad= 6673,
          flux=   21.0,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h21m47.16', '-40.15.47.3'),
          rad=17875,
          flux=   10.5,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h21m47.64', '-43.25.35.2'),
          rad= 7751,
          flux=   41.2,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h21m51.72', '-41.03.15.0'),
          rad=15165,
          flux=   12.0,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h21m51.96', '-45.26.10.3'),
          rad= 5390,
          flux=   15.2,
          bmaj=  25.2,
          bmin= 0.0,
          bpa=  48.9)

addsource(t, dir=dm.direction('J2000', '12h21m55.45', '-43.23.22.0'),
          rad= 7793,
          flux=   12.4,
          bmaj=  28.4,
          bmin= 0.0,
          bpa=  71.7)

addsource(t, dir=dm.direction('J2000', '12h21m55.93', '-45.47.36.6'),
          rad= 5843,
          flux=  101.7,
          bmaj=  40.0,
          bmin=  16.5,
          bpa=  60.6)

addsource(t, dir=dm.direction('J2000', '12h22m01.39', '-44.06.20.7'),
          rad= 6043,
          flux=   14.5,
          bmaj=  22.9,
          bmin= 0.0,
          bpa= -25.4)

addsource(t, dir=dm.direction('J2000', '12h22m04.11', '-47.14.18.5'),
          rad= 9452,
          flux=   31.1,
          bmaj=  17.6,
          bmin= 0.0,
          bpa=  55.3)

addsource(t, dir=dm.direction('J2000', '12h22m07.40', '-40.13.31.7'),
          rad=17938,
          flux=   88.7,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h22m07.74', '-46.45.18.6'),
          rad= 8013,
          flux=  146.4,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h22m14.61', '-40.28.00.4'),
          rad=17085,
          flux=   30.5,
          bmaj=  18.7,
          bmin= 0.0,
          bpa= -20.5)

addsource(t, dir=dm.direction('J2000', '12h22m14.77', '-43.35.51.8'),
          rad= 7100,
          flux=  142.9,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h22m15.05', '-43.30.08.4'),
          rad= 7348,
          flux=   11.7,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h22m17.05', '-42.43.50.1'),
          rad= 9578,
          flux=   14.1,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h22m18.41', '-47.28.17.0'),
          rad=10099,
          flux=   12.2,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h22m21.74', '-41.16.17.3'),
          rad=14317,
          flux=   15.3,
          bmaj=  25.7,
          bmin= 0.0,
          bpa=  20.9)

addsource(t, dir=dm.direction('J2000', '12h22m21.75', '-41.29.07.9'),
          rad=13596,
          flux=   23.8,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h22m24.18', '-44.23.00.4'),
          rad= 5343,
          flux=   21.2,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h22m25.23', '-48.07.20.0'),
          rad=12171,
          flux=   18.4,
          bmaj=  23.3,
          bmin= 0.0,
          bpa= -37.2)

addsource(t, dir=dm.direction('J2000', '12h22m25.80', '-49.25.35.2'),
          rad=16574,
          flux=   10.4,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h22m31.76', '-42.59.12.0'),
          rad= 8711,
          flux=   14.4,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h22m35.13', '-40.36.21.2'),
          rad=16540,
          flux=  463.5,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h22m35.94', '-47.06.47.2'),
          rad= 8898,
          flux=   11.2,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h22m39.96', '-46.39.40.4'),
          rad= 7542,
          flux=   55.0,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h22m41.75', '-48.00.54.0'),
          rad=11752,
          flux=   62.1,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h22m48.42', '-46.51.04.1'),
          rad= 8040,
          flux=   12.9,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h22m53.26', '-43.09.13.5'),
          rad= 8080,
          flux=   32.0,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h22m54.68', '-40.35.21.4'),
          rad=16535,
          flux=   25.9,
          bmaj=  32.4,
          bmin=  21.6,
          bpa=  72.5)

addsource(t, dir=dm.direction('J2000', '12h22m56.68', '-47.00.44.1'),
          rad= 8478,
          flux=   47.6,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h23m02.03', '-42.34.56.7'),
          rad= 9805,
          flux=   26.3,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h23m04.60', '-46.19.01.0'),
          rad= 6436,
          flux=   19.6,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h23m08.74', '-44.35.29.3'),
          rad= 4618,
          flux=  122.1,
          bmaj=  14.6,
          bmin= 0.0,
          bpa=  53.2)

addsource(t, dir=dm.direction('J2000', '12h23m11.26', '-42.56.57.0'),
          rad= 8598,
          flux=   17.0,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h23m12.50', '-44.36.03.0'),
          rad= 4569,
          flux=  131.7,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h23m19.63', '-47.55.39.8'),
          rad=11316,
          flux=   22.6,
          bmaj=  18.5,
          bmin= 0.0,
          bpa=   0.3)

addsource(t, dir=dm.direction('J2000', '12h23m22.92', '-40.16.30.5'),
          rad=17542,
          flux=   10.6,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h23m24.18', '-49.37.11.6'),
          rad=17091,
          flux=   13.0,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h23m25.10', '-47.04.43.7'),
          rad= 8536,
          flux=  241.5,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h23m25.92', '-40.11.34.4'),
          rad=17821,
          flux=   17.8,
          bmaj=  21.2,
          bmin= 0.0,
          bpa= -48.0)

addsource(t, dir=dm.direction('J2000', '12h23m27.88', '-46.47.06.6'),
          rad= 7618,
          flux=   11.7,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h23m31.79', '-40.35.11.4'),
          rad=16435,
          flux=   53.0,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h23m32.16', '-48.58.38.1'),
          rad=14844,
          flux=   12.3,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h23m33.06', '-47.01.56.2'),
          rad= 8350,
          flux=   34.5,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h23m36.18', '-42.47.35.7'),
          rad= 8959,
          flux=   37.2,
          bmaj=  27.9,
          bmin= 0.0,
          bpa=  47.3)

addsource(t, dir=dm.direction('J2000', '12h23m37.79', '-41.12.17.3'),
          rad=14277,
          flux=   20.9,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h23m39.25', '-46.11.19.7'),
          rad= 5854,
          flux=  360.5,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h23m42.42', '-42.32.46.5'),
          rad= 9730,
          flux=   12.9,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h23m49.65', '-47.50.57.3'),
          rad=10943,
          flux=   14.6,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h23m50.40', '-45.37.11.6'),
          rad= 4492,
          flux=   10.6,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h23m50.62', '-45.46.27.3'),
          rad= 4786,
          flux=   65.1,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h23m50.62', '-42.00.01.3'),
          rad=11516,
          flux=   42.3,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h23m51.45', '-49.27.46.6'),
          rad=16481,
          flux=   16.7,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h23m53.11', '-46.50.49.5'),
          rad= 7671,
          flux=   12.8,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h23m54.21', '-45.29.31.1'),
          rad= 4249,
          flux=   11.5,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h23m54.69', '-43.35.15.4'),
          rad= 6420,
          flux=   32.9,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h23m58.22', '-40.34.02.5'),
          rad=16429,
          flux=   11.2,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h23m58.26', '-48.12.58.3'),
          rad=12156,
          flux=  369.2,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h23m59.64', '-44.55.30.1'),
          rad= 3834,
          flux=   24.5,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h24m00.38', '-43.02.27.8'),
          rad= 8046,
          flux=   31.8,
          bmaj=  46.1,
          bmin= 0.0,
          bpa= -15.6)

addsource(t, dir=dm.direction('J2000', '12h24m00.61', '-44.54.27.0'),
          rad= 3829,
          flux=   27.1,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h24m04.61', '-44.02.47.8'),
          rad= 5120,
          flux=   18.9,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h24m05.19', '-48.53.35.7'),
          rad=14466,
          flux=   29.8,
          bmaj=  72.5,
          bmin= 0.0,
          bpa=  76.8)

addsource(t, dir=dm.direction('J2000', '12h24m05.26', '-41.03.42.9'),
          rad=14688,
          flux=   43.3,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h24m05.38', '-41.12.49.7'),
          rad=14161,
          flux=   21.8,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h24m05.63', '-43.53.29.5'),
          rad= 5506,
          flux=   13.6,
          bmaj=  24.0,
          bmin= 0.0,
          bpa= -45.7)

addsource(t, dir=dm.direction('J2000', '12h24m07.40', '-47.23.43.8'),
          rad= 9365,
          flux=   15.0,
          bmaj=  26.4,
          bmin= 0.0,
          bpa= -10.4)

addsource(t, dir=dm.direction('J2000', '12h24m07.52', '-46.18.05.9'),
          rad= 5967,
          flux=   61.5,
          bmaj=  22.3,
          bmin= 0.0,
          bpa=   4.2)

addsource(t, dir=dm.direction('J2000', '12h24m09.19', '-47.20.26.7'),
          rad= 9177,
          flux=   26.7,
          bmaj=  22.5,
          bmin= 0.0,
          bpa= -57.2)

addsource(t, dir=dm.direction('J2000', '12h24m12.59', '-46.20.15.2'),
          rad= 6036,
          flux=  327.8,
          bmaj=  34.5,
          bmin= 0.0,
          bpa= -66.3)

addsource(t, dir=dm.direction('J2000', '12h24m17.88', '-40.17.27.1'),
          rad=17347,
          flux=   14.0,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h24m19.25', '-42.40.56.6'),
          rad= 9118,
          flux=   12.8,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h24m19.97', '-42.10.36.5'),
          rad=10809,
          flux=   11.0,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h24m24.15', '-40.15.10.9'),
          rad=17465,
          flux=   50.1,
          bmaj=  16.7,
          bmin= 0.0,
          bpa= -65.2)

addsource(t, dir=dm.direction('J2000', '12h24m24.34', '-48.03.27.5'),
          rad=11533,
          flux=   34.8,
          bmaj=  19.3,
          bmin= 0.0,
          bpa=  77.3)

addsource(t, dir=dm.direction('J2000', '12h24m27.29', '-42.06.10.3'),
          rad=11033,
          flux=  618.8,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h24m29.25', '-42.06.57.2'),
          rad=10982,
          flux=  733.5,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h24m29.63', '-41.24.13.2'),
          rad=13431,
          flux=  113.1,
          bmaj=  19.4,
          bmin= 0.0,
          bpa= -46.5)

addsource(t, dir=dm.direction('J2000', '12h24m29.81', '-46.08.17.7'),
          rad= 5367,
          flux=   19.0,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h24m30.85', '-42.25.57.9'),
          rad= 9903,
          flux=  180.2,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h24m32.18', '-40.29.56.5'),
          rad=16582,
          flux=   71.4,
          bmaj=  34.9,
          bmin= 0.0,
          bpa= -32.8)

addsource(t, dir=dm.direction('J2000', '12h24m34.15', '-49.14.18.6'),
          rad=15601,
          flux=   29.5,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h24m34.45', '-47.12.21.8'),
          rad= 8631,
          flux=   13.4,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h24m35.46', '-44.33.11.8'),
          rad= 3811,
          flux=   41.9,
          bmaj=  18.2,
          bmin= 0.0,
          bpa=  47.8)

addsource(t, dir=dm.direction('J2000', '12h24m35.83', '-44.11.09.3'),
          rad= 4536,
          flux=   29.4,
          bmaj=  20.5,
          bmin= 0.0,
          bpa= -36.5)

addsource(t, dir=dm.direction('J2000', '12h24m35.96', '-49.35.31.3'),
          rad=16837,
          flux=   15.0,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h24m40.89', '-45.55.20.8'),
          rad= 4722,
          flux=   63.3,
          bmaj=  16.0,
          bmin= 0.0,
          bpa= -36.7)

addsource(t, dir=dm.direction('J2000', '12h24m45.58', '-43.18.44.2'),
          rad= 6953,
          flux=   11.6,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h24m50.25', '-49.31.41.6'),
          rad=16585,
          flux=   29.6,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h24m51.52', '-43.35.20.8'),
          rad= 6063,
          flux=  366.0,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h24m54.47', '-40.21.20.8'),
          rad=17035,
          flux=   33.4,
          bmaj=  50.0,
          bmin= 0.0,
          bpa= -67.9)

addsource(t, dir=dm.direction('J2000', '12h24m54.81', '-40.54.29.5'),
          rad=15093,
          flux=   57.3,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h24m54.86', '-43.13.11.9'),
          rad= 7200,
          flux=   23.4,
          bmaj=  23.6,
          bmin= 0.0,
          bpa= -73.2)

addsource(t, dir=dm.direction('J2000', '12h24m57.14', '-49.09.08.9'),
          rad=15251,
          flux=  101.3,
          bmaj=  25.4,
          bmin= 0.0,
          bpa= -24.1)

addsource(t, dir=dm.direction('J2000', '12h24m58.21', '-47.37.14.0'),
          rad= 9934,
          flux=   22.1,
          bmaj=  23.7,
          bmin= 0.0,
          bpa=  -5.8)

addsource(t, dir=dm.direction('J2000', '12h24m59.96', '-40.48.35.6'),
          rad=15425,
          flux=  168.2,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h25m00.18', '-40.09.27.7'),
          rad=17722,
          flux=   12.9,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h25m03.94', '-43.12.13.7'),
          rad= 7208,
          flux=   21.2,
          bmaj=  45.9,
          bmin= 0.0,
          bpa=  85.1)

addsource(t, dir=dm.direction('J2000', '12h25m06.61', '-48.34.36.3'),
          rad=13215,
          flux=   76.1,
          bmaj=  14.5,
          bmin= 0.0,
          bpa=  40.8)

addsource(t, dir=dm.direction('J2000', '12h25m11.85', '-41.34.06.7'),
          rad=12739,
          flux=   26.5,
          bmaj=  24.5,
          bmin= 0.0,
          bpa=  69.3)

addsource(t, dir=dm.direction('J2000', '12h25m12.07', '-41.27.21.9'),
          rad=13131,
          flux=   12.8,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h25m12.39', '-45.23.29.2'),
          rad= 3351,
          flux=   96.2,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h25m12.51', '-45.42.14.3'),
          rad= 3950,
          flux=   23.4,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h25m13.39', '-42.42.33.5'),
          rad= 8807,
          flux=   15.4,
          bmaj=  32.6,
          bmin= 0.0,
          bpa= -74.9)

addsource(t, dir=dm.direction('J2000', '12h25m14.39', '-43.13.40.0'),
          rad= 7081,
          flux=   14.2,
          bmaj=  24.8,
          bmin= 0.0,
          bpa= -10.2)

addsource(t, dir=dm.direction('J2000', '12h25m16.56', '-46.22.12.2'),
          rad= 5757,
          flux=   18.4,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h25m17.89', '-45.41.46.5'),
          rad= 3889,
          flux=   20.4,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h25m21.83', '-45.58.20.5'),
          rad= 4561,
          flux=   20.4,
          bmaj=  66.6,
          bmin=  24.9,
          bpa= -59.7)

addsource(t, dir=dm.direction('J2000', '12h25m24.68', '-45.15.03.6'),
          rad= 3051,
          flux=  346.7,
          bmaj=  43.7,
          bmin= 0.0,
          bpa=  10.7)

addsource(t, dir=dm.direction('J2000', '12h25m24.76', '-43.34.39.2'),
          rad= 5911,
          flux=   25.7,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h25m25.45', '-40.27.39.2'),
          rad=16600,
          flux=   13.4,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h25m28.30', '-41.26.07.0'),
          rad=13163,
          flux=   76.6,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h25m29.45', '-45.31.32.4'),
          rad= 3426,
          flux=   59.0,
          bmaj=  21.0,
          bmin= 0.0,
          bpa=  31.2)

addsource(t, dir=dm.direction('J2000', '12h25m29.69', '-40.35.48.1'),
          rad=16112,
          flux=   13.5,
          bmaj=  31.9,
          bmin= 0.0,
          bpa=  73.0)

addsource(t, dir=dm.direction('J2000', '12h25m30.37', '-46.29.51.7'),
          rad= 6085,
          flux=   21.8,
          bmaj=  27.8,
          bmin= 0.0,
          bpa=   4.1)

addsource(t, dir=dm.direction('J2000', '12h25m34.51', '-40.46.14.3'),
          rad=15488,
          flux=   16.0,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h25m35.45', '-45.41.04.6'),
          rad= 3722,
          flux=   19.5,
          bmaj=  50.4,
          bmin=  30.2,
          bpa= -69.3)

addsource(t, dir=dm.direction('J2000', '12h25m36.37', '-41.05.09.0'),
          rad=14372,
          flux=   19.8,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h25m38.10', '-44.48.36.8'),
          rad= 2865,
          flux=   29.8,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h25m38.84', '-40.18.05.8'),
          rad=17138,
          flux=   17.2,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h25m39.08', '-42.59.00.2'),
          rad= 7785,
          flux=   26.0,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h25m41.37', '-49.04.34.2'),
          rad=14897,
          flux=   44.8,
          bmaj=  63.6,
          bmin=  22.6,
          bpa=  82.6)

addsource(t, dir=dm.direction('J2000', '12h25m41.66', '-49.48.46.5'),
          rad=17502,
          flux=   24.9,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h25m44.84', '-47.04.01.9'),
          rad= 7900,
          flux=   13.3,
          bmaj=  36.0,
          bmin= 0.0,
          bpa= -40.1)

addsource(t, dir=dm.direction('J2000', '12h25m46.77', '-41.40.56.7'),
          rad=12251,
          flux=   50.3,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h25m48.79', '-43.35.32.9'),
          rad= 5739,
          flux=   12.0,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h25m52.59', '-42.27.49.5'),
          rad= 9513,
          flux=   21.0,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h25m53.70', '-48.12.28.0'),
          rad=11817,
          flux=   10.6,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h25m54.78', '-44.56.22.4'),
          rad= 2611,
          flux=  429.3,
          bmaj=  26.0,
          bmin= 0.0,
          bpa= -65.1)

addsource(t, dir=dm.direction('J2000', '12h25m55.06', '-42.36.50.9'),
          rad= 8986,
          flux=  115.6,
          bmaj=  17.2,
          bmin= 0.0,
          bpa= -41.8)

addsource(t, dir=dm.direction('J2000', '12h25m55.61', '-47.16.21.2'),
          rad= 8564,
          flux=   12.5,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h25m56.09', '-45.38.51.6'),
          rad= 3472,
          flux=   15.7,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h25m58.06', '-43.33.12.1'),
          rad= 5819,
          flux=   12.0,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h25m59.91', '-47.09.30.0'),
          rad= 8159,
          flux=   11.1,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h26m00.05', '-40.51.17.7'),
          rad=15139,
          flux=   19.7,
          bmaj=  24.8,
          bmin= 0.0,
          bpa=   0.8)

addsource(t, dir=dm.direction('J2000', '12h26m00.55', '-45.51.53.5'),
          rad= 4005,
          flux=   79.9,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h26m07.92', '-47.37.00.3'),
          rad= 9719,
          flux=  160.4,
          bmaj=  23.8,
          bmin= 0.0,
          bpa=  36.4)

addsource(t, dir=dm.direction('J2000', '12h26m08.33', '-44.56.43.5'),
          rad= 2466,
          flux=   11.2,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h26m09.25', '-45.58.50.2'),
          rad= 4283,
          flux=   27.1,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h26m11.92', '-40.47.33.7'),
          rad=15338,
          flux=   14.1,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h26m13.12', '-47.18.17.6'),
          rad= 8623,
          flux=   16.1,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h26m13.17', '-42.28.39.0'),
          rad= 9404,
          flux=   17.5,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h26m22.55', '-48.36.40.4'),
          rad=13181,
          flux=   16.4,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h26m24.62', '-46.53.28.8'),
          rad= 7168,
          flux=  117.3,
          bmaj=  32.6,
          bmin= 0.0,
          bpa=   7.4)

addsource(t, dir=dm.direction('J2000', '12h26m27.78', '-48.12.48.2'),
          rad=11767,
          flux=   13.6,
          bmaj=  35.1,
          bmin= 0.0,
          bpa= -89.7)

addsource(t, dir=dm.direction('J2000', '12h26m30.04', '-49.12.08.0'),
          rad=15265,
          flux=   23.4,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h26m34.73', '-40.27.01.5'),
          rad=16516,
          flux=   10.9,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h26m38.62', '-42.32.39.8'),
          rad= 9102,
          flux=  238.2,
          bmaj=  39.9,
          bmin= 0.0,
          bpa= -89.0)

addsource(t, dir=dm.direction('J2000', '12h26m39.92', '-43.05.14.5'),
          rad= 7214,
          flux=   35.5,
          bmaj=  72.1,
          bmin=  18.4,
          bpa= -64.4)

addsource(t, dir=dm.direction('J2000', '12h26m41.05', '-48.47.23.5'),
          rad=13785,
          flux=   59.6,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h26m42.95', '-47.50.30.7'),
          rad=10427,
          flux=   12.6,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h26m43.37', '-44.51.17.2'),
          rad= 2153,
          flux=   10.7,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h26m43.44', '-46.16.38.9'),
          rad= 5039,
          flux=   20.0,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h26m43.67', '-45.32.01.9'),
          rad= 2826,
          flux=   14.2,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h26m44.79', '-40.35.48.6'),
          rad=15980,
          flux=   19.0,
          bmaj=  45.7,
          bmin= 0.0,
          bpa=  67.0)

addsource(t, dir=dm.direction('J2000', '12h26m46.84', '-48.31.29.4'),
          rad=12835,
          flux=   41.4,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h26m48.54', '-43.32.19.6'),
          rad= 5647,
          flux=   15.8,
          bmaj=  37.3,
          bmin= 0.0,
          bpa= -29.9)

addsource(t, dir=dm.direction('J2000', '12h26m51.82', '-43.43.20.7'),
          rad= 5022,
          flux=   42.5,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h26m51.98', '-41.00.26.8'),
          rad=14508,
          flux=   27.1,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h26m52.51', '-42.18.19.0'),
          rad= 9908,
          flux=   65.5,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h26m52.74', '-49.22.40.6'),
          rad=15860,
          flux=   17.6,
          bmaj=  33.8,
          bmin=  22.1,
          bpa= -74.2)

addsource(t, dir=dm.direction('J2000', '12h26m53.88', '-40.10.44.8'),
          rad=17455,
          flux=   26.4,
          bmaj=  83.5,
          bmin= 0.0,
          bpa=  36.5)

addsource(t, dir=dm.direction('J2000', '12h26m54.34', '-43.22.45.5'),
          rad= 6166,
          flux=   35.7,
          bmaj=  16.0,
          bmin= 0.0,
          bpa=  37.3)

addsource(t, dir=dm.direction('J2000', '12h26m55.67', '-44.05.07.6'),
          rad= 3837,
          flux=   23.1,
          bmaj=  34.1,
          bmin= 0.0,
          bpa=  76.9)

addsource(t, dir=dm.direction('J2000', '12h26m57.01', '-48.18.35.3'),
          rad=12056,
          flux=   12.8,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h26m57.94', '-43.40.58.6'),
          rad= 5127,
          flux=  180.6,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h27m05.79', '-40.30.21.4'),
          rad=16275,
          flux=   17.1,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h27m10.36', '-44.07.41.2'),
          rad= 3625,
          flux=   36.4,
          bmaj=  20.7,
          bmin= 0.0,
          bpa= -17.8)

addsource(t, dir=dm.direction('J2000', '12h27m10.98', '-46.48.44.9'),
          rad= 6758,
          flux=   14.2,
          bmaj=  39.2,
          bmin= 0.0,
          bpa= -86.6)

addsource(t, dir=dm.direction('J2000', '12h27m15.45', '-48.11.07.5'),
          rad=11586,
          flux=   71.9,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h27m19.31', '-43.58.11.2'),
          rad= 4088,
          flux=   27.0,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h27m21.11', '-42.59.23.7'),
          rad= 7435,
          flux=   44.7,
          bmaj=  34.2,
          bmin= 0.0,
          bpa= -78.8)

addsource(t, dir=dm.direction('J2000', '12h27m22.87', '-46.32.11.1'),
          rad= 5769,
          flux=   18.1,
          bmaj=  19.6,
          bmin= 0.0,
          bpa=  83.7)

addsource(t, dir=dm.direction('J2000', '12h27m26.06', '-40.29.11.2'),
          rad=16320,
          flux=  104.9,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h27m26.22', '-42.59.22.1'),
          rad= 7424,
          flux=   37.0,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h27m27.31', '-43.00.54.6'),
          rad= 7331,
          flux=   14.2,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h27m31.00', '-49.28.44.5'),
          rad=16179,
          flux=  121.3,
          bmaj=  25.7,
          bmin= 0.0,
          bpa=  -5.9)

addsource(t, dir=dm.direction('J2000', '12h27m34.81', '-44.58.32.0'),
          rad= 1543,
          flux=   96.6,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h27m38.11', '-40.09.38.0'),
          rad=17471,
          flux=   18.9,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h27m42.04', '-44.00.42.1'),
          rad= 3852,
          flux=  383.7,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h27m43.18', '-45.04.59.6'),
          rad= 1481,
          flux=   21.4,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h27m45.03', '-49.14.33.0'),
          rad=15321,
          flux=   11.2,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h27m52.59', '-49.40.36.7'),
          rad=16867,
          flux=   22.2,
          bmaj=  79.0,
          bmin= 0.0,
          bpa=  60.8)

addsource(t, dir=dm.direction('J2000', '12h27m54.16', '-40.57.13.5'),
          rad=14619,
          flux=   21.5,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h27m55.11', '-48.42.34.7'),
          rad=13406,
          flux=   32.6,
          bmaj=  28.8,
          bmin= 0.0,
          bpa= -89.8)

addsource(t, dir=dm.direction('J2000', '12h27m55.68', '-49.32.56.3'),
          rad=16408,
          flux=  235.2,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h27m55.78', '-43.59.56.6'),
          rad= 3840,
          flux=   12.3,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h27m56.38', '-40.59.01.4'),
          rad=14510,
          flux=  118.7,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h27m56.91', '-45.54.57.5'),
          rad= 3542,
          flux=   15.6,
          bmaj=  30.2,
          bmin= 0.0,
          bpa=  34.3)

addsource(t, dir=dm.direction('J2000', '12h27m58.52', '-40.56.04.5'),
          rad=14684,
          flux=   20.4,
          bmaj=  19.2,
          bmin= 0.0,
          bpa= -56.6)

addsource(t, dir=dm.direction('J2000', '12h28m01.12', '-49.30.52.8'),
          rad=16281,
          flux=   13.2,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h28m03.35', '-46.51.51.1'),
          rad= 6819,
          flux=   21.3,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h28m05.59', '-48.19.28.9'),
          rad=12020,
          flux=   16.3,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h28m09.88', '-43.45.17.4'),
          rad= 4635,
          flux=   12.5,
          bmaj=  27.1,
          bmin= 0.0,
          bpa=  17.2)

addsource(t, dir=dm.direction('J2000', '12h28m10.72', '-44.46.04.8'),
          rad= 1431,
          flux=   14.9,
          bmaj=  31.8,
          bmin= 0.0,
          bpa= -56.6)

addsource(t, dir=dm.direction('J2000', '12h28m11.06', '-44.05.34.1'),
          rad= 3467,
          flux=   48.7,
          bmaj=  44.9,
          bmin=  29.3,
          bpa=  36.0)

addsource(t, dir=dm.direction('J2000', '12h28m12.17', '-44.06.57.0'),
          rad= 3385,
          flux=   29.6,
          bmaj=  54.8,
          bmin= 0.0,
          bpa=  61.4)

addsource(t, dir=dm.direction('J2000', '12h28m12.21', '-49.48.22.5'),
          rad=17317,
          flux=   51.0,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h28m13.24', '-40.04.13.8'),
          rad=17763,
          flux=  124.4,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h28m14.66', '-44.06.33.4'),
          rad= 3398,
          flux=  112.4,
          bmaj= 127.0,
          bmin=  83.4,
          bpa= -54.4)

addsource(t, dir=dm.direction('J2000', '12h28m16.69', '-46.20.09.7'),
          rad= 4930,
          flux=   17.3,
          bmaj=  20.5,
          bmin= 0.0,
          bpa=  79.5)

addsource(t, dir=dm.direction('J2000', '12h28m18.30', '-44.59.38.5'),
          rad= 1079,
          flux=   15.0,
          bmaj=  36.5,
          bmin= 0.0,
          bpa= -15.0)

addsource(t, dir=dm.direction('J2000', '12h28m19.68', '-41.02.33.0'),
          rad=14278,
          flux=   13.7,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h28m21.29', '-41.03.48.5'),
          rad=14201,
          flux=   47.4,
          bmaj=  17.4,
          bmin= 0.0,
          bpa= -30.8)

addsource(t, dir=dm.direction('J2000', '12h28m29.98', '-42.16.32.6'),
          rad= 9852,
          flux=   13.6,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h28m37.21', '-46.02.28.8'),
          rad= 3848,
          flux=   10.7,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h28m37.69', '-40.57.55.1'),
          rad=14541,
          flux=   13.0,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h28m42.15', '-43.39.33.9'),
          rad= 4897,
          flux=   18.5,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h28m47.95', '-44.47.14.5'),
          rad= 1083,
          flux=   19.6,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h28m49.70', '-43.42.23.6'),
          rad= 4717,
          flux=   10.7,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h28m49.72', '-49.56.54.9'),
          rad=17807,
          flux=   52.4,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h28m50.20', '-44.47.19.5'),
          rad= 1062,
          flux=  294.6,
          bmaj=  63.2,
          bmin=  18.3,
          bpa=  72.8)

addsource(t, dir=dm.direction('J2000', '12h28m50.84', '-40.06.34.2'),
          rad=17601,
          flux=   13.5,
          bmaj=  31.8,
          bmin= 0.0,
          bpa= -82.0)

addsource(t, dir=dm.direction('J2000', '12h28m51.78', '-48.58.01.7'),
          rad=14287,
          flux=  387.0,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h28m51.99', '-47.21.32.3'),
          rad= 8519,
          flux=   25.5,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h28m52.47', '-49.04.31.9'),
          rad=14676,
          flux=  105.2,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h28m55.08', '-41.30.32.2'),
          rad=12580,
          flux=   38.3,
          bmaj=  21.5,
          bmin=  17.1,
          bpa= -20.8)

addsource(t, dir=dm.direction('J2000', '12h28m56.68', '-49.02.15.6'),
          rad=14538,
          flux=   31.3,
          bmaj=  17.2,
          bmin= 0.0,
          bpa=  58.4)

addsource(t, dir=dm.direction('J2000', '12h28m56.99', '-44.20.29.3'),
          rad= 2464,
          flux=   54.2,
          bmaj=  63.3,
          bmin=  18.1,
          bpa=  18.0)

addsource(t, dir=dm.direction('J2000', '12h28m57.36', '-42.50.05.5'),
          rad= 7822,
          flux=   27.6,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h29m01.38', '-47.36.55.6'),
          rad= 9432,
          flux=   18.8,
          bmaj=  47.6,
          bmin= 0.0,
          bpa=  83.1)

addsource(t, dir=dm.direction('J2000', '12h29m01.83', '-43.34.09.1'),
          rad= 5188,
          flux=   11.9,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h29m03.28', '-47.07.19.6'),
          rad= 7661,
          flux=   20.3,
          bmaj=  23.6,
          bmin= 0.0,
          bpa=  20.7)

addsource(t, dir=dm.direction('J2000', '12h29m04.89', '-43.30.28.1'),
          rad= 5404,
          flux=   14.7,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h29m08.32', '-42.58.59.2'),
          rad= 7281,
          flux=   13.4,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h29m10.84', '-41.13.03.7'),
          rad=13617,
          flux=   18.9,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h29m20.80', '-44.28.51.1'),
          rad= 1915,
          flux=   30.6,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h29m20.81', '-41.38.35.3'),
          rad=12085,
          flux=   13.0,
          bmaj=  27.1,
          bmin= 0.0,
          bpa=  48.1)

addsource(t, dir=dm.direction('J2000', '12h29m21.30', '-49.07.33.0'),
          rad=14845,
          flux=   15.7,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h29m22.83', '-40.30.15.3'),
          rad=16173,
          flux=   17.0,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h29m23.03', '-46.23.56.3'),
          rad= 5051,
          flux=   24.6,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h29m26.57', '-42.18.55.3'),
          rad= 9668,
          flux=   17.6,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h29m28.50', '-41.15.19.6'),
          rad=13475,
          flux=   13.1,
          bmaj=  26.0,
          bmin= 0.0,
          bpa= -58.0)

addsource(t, dir=dm.direction('J2000', '12h29m36.70', '-45.33.08.8'),
          rad= 2004,
          flux=   13.9,
          bmaj=  25.1,
          bmin= 0.0,
          bpa=  -9.6)

addsource(t, dir=dm.direction('J2000', '12h29m39.89', '-43.53.19.3'),
          rad= 4006,
          flux=   33.6,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h29m44.37', '-42.05.04.9'),
          rad=10492,
          flux=   13.9,
          bmaj=  25.2,
          bmin= 0.0,
          bpa=  43.4)

addsource(t, dir=dm.direction('J2000', '12h29m45.54', '-43.28.17.3'),
          rad= 5504,
          flux=   53.6,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h29m46.29', '-48.53.59.6'),
          rad=14029,
          flux=   43.8,
          bmaj=  16.2,
          bmin= 0.0,
          bpa=  68.4)

addsource(t, dir=dm.direction('J2000', '12h29m51.27', '-46.01.12.5'),
          rad= 3673,
          flux=   16.1,
          bmaj=  45.7,
          bmin= 0.0,
          bpa=  28.4)

addsource(t, dir=dm.direction('J2000', '12h29m56.74', '-45.23.29.2'),
          rad= 1410,
          flux=   18.4,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h29m57.56', '-40.21.28.0'),
          rad=16694,
          flux=   13.6,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h29m59.25', '-45.04.26.6'),
          rad=  267,
          flux=   11.9,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h30m00.81', '-44.36.58.9'),
          rad= 1381,
          flux=   11.5,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h30m03.78', '-45.19.24.2'),
          rad= 1165,
          flux=   65.6,
          bmaj=  17.3,
          bmin= 0.0,
          bpa= -41.0)

addsource(t, dir=dm.direction('J2000', '12h30m07.30', '-44.18.16.8'),
          rad= 2504,
          flux=   15.2,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h30m08.46', '-42.33.12.3'),
          rad= 8805,
          flux=   15.5,
          bmaj=  43.1,
          bmin= 0.0,
          bpa=  80.1)

addsource(t, dir=dm.direction('J2000', '12h30m11.80', '-47.00.23.0'),
          rad= 7223,
          flux=   93.6,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h30m15.56', '-49.40.49.7'),
          rad=16832,
          flux=   30.7,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h30m16.94', '-43.49.17.3'),
          rad= 4246,
          flux=   44.2,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h30m17.63', '-46.21.15.7'),
          rad= 4879,
          flux=   12.2,
          bmaj=  29.8,
          bmin= 0.0,
          bpa= -22.6)

addsource(t, dir=dm.direction('J2000', '12h30m18.38', '-46.06.00.6'),
          rad= 3965,
          flux=   14.3,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h30m19.34', '-44.28.40.2'),
          rad= 1891,
          flux=  111.5,
          bmaj=  17.4,
          bmin=  14.6,
          bpa=  35.6)

addsource(t, dir=dm.direction('J2000', '12h30m23.56', '-49.09.26.6'),
          rad=14955,
          flux=   11.1,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h30m29.94', '-44.11.14.7'),
          rad= 2943,
          flux=   13.8,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h30m30.07', '-41.40.38.8'),
          rad=11959,
          flux=  114.9,
          bmaj= 127.1,
          bmin=  99.4,
          bpa=   1.2)

addsource(t, dir=dm.direction('J2000', '12h30m32.25', '-47.28.34.9'),
          rad= 8918,
          flux=   53.4,
          bmaj=  14.9,
          bmin= 0.0,
          bpa=  51.0)

addsource(t, dir=dm.direction('J2000', '12h30m32.65', '-49.25.08.9'),
          rad=15897,
          flux=   30.6,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h30m33.00', '-49.10.15.8'),
          rad=15006,
          flux=   13.1,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h30m35.01', '-41.38.51.4'),
          rad=12068,
          flux=  648.6,
          bmaj= 127.2,
          bmin=  67.7,
          bpa= -59.8)

addsource(t, dir=dm.direction('J2000', '12h30m37.76', '-48.51.40.0'),
          rad=13895,
          flux=   31.2,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h30m43.68', '-43.26.09.0'),
          rad= 5650,
          flux=   12.1,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h30m44.87', '-41.37.57.8'),
          rad=12125,
          flux=   44.4,
          bmaj=  99.4,
          bmin=  44.3,
          bpa=  25.5)

addsource(t, dir=dm.direction('J2000', '12h30m48.79', '-48.57.58.6'),
          rad=14276,
          flux=  196.0,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h30m53.71', '-45.53.27.8'),
          rad= 3257,
          flux=   17.5,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h30m55.35', '-41.47.05.2'),
          rad=11584,
          flux=   14.1,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h30m55.98', '-40.26.39.7'),
          rad=16395,
          flux=   13.4,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h30m59.32', '-45.38.35.8'),
          rad= 2399,
          flux=   11.2,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h30m59.61', '-41.20.51.9'),
          rad=13155,
          flux=   13.5,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h31m01.84', '-47.56.53.3'),
          rad=10628,
          flux=   44.5,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h31m02.95', '-49.34.39.0'),
          rad=16474,
          flux=   10.9,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h31m03.69', '-45.26.27.4'),
          rad= 1724,
          flux=   25.5,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h31m06.56', '-41.09.59.0'),
          rad=13810,
          flux=  180.1,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h31m07.18', '-47.59.07.3'),
          rad=10765,
          flux=  144.8,
          bmaj=  27.5,
          bmin= 0.0,
          bpa= -36.4)

addsource(t, dir=dm.direction('J2000', '12h31m12.23', '-42.06.10.2'),
          rad=10455,
          flux=   62.7,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h31m15.30', '-47.47.26.4'),
          rad=10072,
          flux=   14.4,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h31m16.46', '-48.57.09.9'),
          rad=14240,
          flux=   17.9,
          bmaj=  30.7,
          bmin= 0.0,
          bpa= -48.1)

addsource(t, dir=dm.direction('J2000', '12h31m17.21', '-49.03.21.5'),
          rad=14611,
          flux=   15.1,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h31m18.96', '-49.22.25.5'),
          rad=15751,
          flux=   13.2,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h31m21.20', '-45.28.28.1'),
          rad= 1911,
          flux=   37.8,
          bmaj=  15.8,
          bmin= 0.0,
          bpa= -49.8)

addsource(t, dir=dm.direction('J2000', '12h31m22.48', '-41.46.33.4'),
          rad=11635,
          flux=   26.7,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h31m24.57', '-45.11.20.7'),
          rad= 1125,
          flux=   45.3,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h31m25.04', '-46.02.47.9'),
          rad= 3872,
          flux=   27.6,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h31m28.58', '-41.49.23.0'),
          rad=11472,
          flux=  104.0,
          bmaj=  21.1,
          bmin= 0.0,
          bpa=  45.7)

addsource(t, dir=dm.direction('J2000', '12h31m29.27', '-49.15.37.9'),
          rad=15351,
          flux=   57.1,
          bmaj=  51.8,
          bmin=  21.0,
          bpa= -14.9)

addsource(t, dir=dm.direction('J2000', '12h31m33.28', '-43.59.50.2'),
          rad= 3745,
          flux=   25.6,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h31m35.10', '-48.02.13.9'),
          rad=10973,
          flux=   20.7,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h31m35.79', '-41.41.14.4'),
          rad=11965,
          flux=   22.7,
          bmaj=  19.7,
          bmin= 0.0,
          bpa= -10.9)

addsource(t, dir=dm.direction('J2000', '12h31m40.33', '-40.17.31.5'),
          rad=16965,
          flux=  353.4,
          bmaj=  15.2,
          bmin= 0.0,
          bpa=  88.8)

addsource(t, dir=dm.direction('J2000', '12h31m42.53', '-44.31.06.1'),
          rad= 2049,
          flux=   13.7,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h31m43.00', '-47.59.20.1'),
          rad=10808,
          flux=   15.3,
          bmaj=  34.3,
          bmin= 0.0,
          bpa= -58.6)

addsource(t, dir=dm.direction('J2000', '12h31m43.80', '-41.37.37.7'),
          rad=12188,
          flux=   39.1,
          bmaj=  26.9,
          bmin= 0.0,
          bpa= -59.1)

addsource(t, dir=dm.direction('J2000', '12h31m43.96', '-46.56.48.8'),
          rad= 7091,
          flux=   27.7,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h31m44.59', '-42.09.56.3'),
          rad=10262,
          flux=   19.2,
          bmaj=  23.9,
          bmin= 0.0,
          bpa= -38.7)

addsource(t, dir=dm.direction('J2000', '12h31m47.90', '-42.34.08.9'),
          rad= 8826,
          flux=   26.7,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h31m49.64', '-40.21.33.1'),
          rad=16732,
          flux=   16.5,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h31m51.99', '-44.49.40.2'),
          rad= 1341,
          flux=   17.7,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h31m55.28', '-45.54.37.6'),
          rad= 3495,
          flux=   23.4,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h31m57.40', '-45.59.49.2'),
          rad= 3795,
          flux=   13.8,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h32m05.86', '-49.10.26.7'),
          rad=15068,
          flux=   11.1,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h32m09.13', '-43.09.34.6'),
          rad= 6769,
          flux=   16.5,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h32m12.88', '-48.41.43.6'),
          rad=13364,
          flux=   17.2,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h32m12.89', '-41.29.33.6'),
          rad=12702,
          flux=   27.1,
          bmaj=  26.6,
          bmin= 0.0,
          bpa= -39.1)

addsource(t, dir=dm.direction('J2000', '12h32m16.51', '-41.50.57.6'),
          rad=11434,
          flux=   13.4,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h32m22.53', '-41.34.34.6'),
          rad=12416,
          flux=   37.9,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h32m23.04', '-45.18.06.4'),
          rad= 1863,
          flux=   11.7,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h32m24.39', '-45.24.27.3'),
          rad= 2117,
          flux=   18.9,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h32m28.63', '-43.07.09.9'),
          rad= 6956,
          flux=   12.6,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h32m29.33', '-46.30.16.6'),
          rad= 5637,
          flux=   39.7,
          bmaj=  18.4,
          bmin= 0.0,
          bpa= -20.8)

addsource(t, dir=dm.direction('J2000', '12h32m30.35', '-44.07.51.9'),
          rad= 3516,
          flux=   37.3,
          bmaj=  29.6,
          bmin=  21.4,
          bpa=  29.4)

addsource(t, dir=dm.direction('J2000', '12h32m31.00', '-44.19.47.5'),
          rad= 2901,
          flux=   40.4,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h32m32.52', '-44.09.06.5'),
          rad= 3461,
          flux=  131.4,
          bmaj=  30.1,
          bmin=  17.6,
          bpa=  14.6)

addsource(t, dir=dm.direction('J2000', '12h32m34.67', '-45.46.05.9'),
          rad= 3210,
          flux=   13.7,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h32m34.78', '-48.21.33.4'),
          rad=12191,
          flux=  325.9,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h32m37.63', '-48.17.05.8'),
          rad=11930,
          flux=  189.6,
          bmaj=  18.8,
          bmin= 0.0,
          bpa= -84.6)

addsource(t, dir=dm.direction('J2000', '12h32m38.41', '-47.53.44.2'),
          rad=10547,
          flux=   32.3,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h32m39.62', '-49.49.24.3'),
          rad=17419,
          flux=   43.7,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h32m39.86', '-46.21.43.6'),
          rad= 5181,
          flux=   18.5,
          bmaj=  36.2,
          bmin=  25.6,
          bpa=   1.7)

addsource(t, dir=dm.direction('J2000', '12h32m42.70', '-41.17.19.9'),
          rad=13469,
          flux=   42.0,
          bmaj=  17.3,
          bmin= 0.0,
          bpa= -23.4)

addsource(t, dir=dm.direction('J2000', '12h32m45.86', '-48.20.56.5'),
          rad=12170,
          flux=   30.6,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h32m49.74', '-43.55.48.5'),
          rad= 4258,
          flux=   22.4,
          bmaj=  37.6,
          bmin= 0.0,
          bpa=   0.4)

addsource(t, dir=dm.direction('J2000', '12h32m50.47', '-42.33.35.2'),
          rad= 8974,
          flux=   45.5,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h32m50.74', '-41.57.40.2'),
          rad=11091,
          flux=   26.8,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h32m51.73', '-45.55.00.6'),
          rad= 3763,
          flux=   71.5,
          bmaj=  54.8,
          bmin= 0.0,
          bpa= -42.4)

addsource(t, dir=dm.direction('J2000', '12h32m52.45', '-40.53.42.3'),
          rad=14885,
          flux=  105.7,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h32m53.74', '-45.42.09.3'),
          rad= 3123,
          flux=   71.8,
          bmaj=  72.1,
          bmin= 0.0,
          bpa= -82.7)

addsource(t, dir=dm.direction('J2000', '12h32m54.35', '-49.36.04.0'),
          rad=16640,
          flux=   18.1,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h32m55.22', '-42.40.16.0'),
          rad= 8593,
          flux=   16.4,
          bmaj=  25.5,
          bmin= 0.0,
          bpa=  45.3)

addsource(t, dir=dm.direction('J2000', '12h32m56.70', '-47.27.44.7'),
          rad= 9049,
          flux=   10.8,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h32m57.75', '-42.48.15.9'),
          rad= 8132,
          flux=   36.5,
          bmaj=  17.0,
          bmin= 0.0,
          bpa= -71.6)

addsource(t, dir=dm.direction('J2000', '12h32m59.19', '-48.08.52.7'),
          rad=11476,
          flux=   62.3,
          bmaj=  16.8,
          bmin= 0.0,
          bpa= -68.2)

addsource(t, dir=dm.direction('J2000', '12h32m59.41', '-49.07.06.2'),
          rad=14926,
          flux=   29.1,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h32m59.62', '-42.22.39.4'),
          rad= 9636,
          flux=   18.4,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h32m59.90', '-40.07.21.3'),
          rad=17649,
          flux=   13.9,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h33m02.54', '-42.12.33.8'),
          rad=10236,
          flux=   33.7,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h33m03.11', '-45.26.46.7'),
          rad= 2515,
          flux=   15.3,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h33m03.31', '-49.01.20.3'),
          rad=14589,
          flux=   44.6,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h33m03.88', '-46.06.53.4'),
          rad= 4454,
          flux=   18.3,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h33m06.65', '-45.29.01.0'),
          rad= 2630,
          flux=  169.9,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h33m08.38', '-47.00.05.7'),
          rad= 7466,
          flux=   30.0,
          bmaj=  39.5,
          bmin= 0.0,
          bpa=  12.6)

addsource(t, dir=dm.direction('J2000', '12h33m10.60', '-46.32.02.1'),
          rad= 5870,
          flux=   16.5,
          bmaj=  28.4,
          bmin= 0.0,
          bpa= -19.9)

addsource(t, dir=dm.direction('J2000', '12h33m14.16', '-43.39.13.0'),
          rad= 5275,
          flux=  212.1,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h33m15.82', '-49.28.35.9'),
          rad=16222,
          flux=   20.8,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h33m16.43', '-48.42.03.4'),
          rad=13465,
          flux=   39.2,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h33m22.87', '-47.38.43.0'),
          rad= 9748,
          flux=   23.4,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h33m32.47', '-43.38.49.1'),
          rad= 5377,
          flux=   27.8,
          bmaj=  21.2,
          bmin= 0.0,
          bpa= -13.5)

addsource(t, dir=dm.direction('J2000', '12h33m33.95', '-40.53.04.6'),
          rad=14987,
          flux=   59.9,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h33m34.05', '-43.33.59.3'),
          rad= 5649,
          flux=   60.0,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h33m34.54', '-49.53.10.5'),
          rad=17702,
          flux=   29.3,
          bmaj=  19.9,
          bmin= 0.0,
          bpa=  51.6)

addsource(t, dir=dm.direction('J2000', '12h33m35.05', '-48.54.22.5'),
          rad=14222,
          flux=  120.0,
          bmaj=  25.2,
          bmin= 0.0,
          bpa= -43.9)

addsource(t, dir=dm.direction('J2000', '12h33m37.81', '-43.12.31.5'),
          rad= 6861,
          flux=   87.8,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h33m38.03', '-41.55.46.4'),
          rad=11300,
          flux=   12.7,
          bmaj=  26.4,
          bmin= 0.0,
          bpa=  43.4)

addsource(t, dir=dm.direction('J2000', '12h33m39.22', '-43.23.26.5'),
          rad= 6254,
          flux=   44.7,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h33m40.24', '-43.01.43.6'),
          rad= 7482,
          flux=   12.2,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h33m42.64', '-47.54.16.4'),
          rad=10702,
          flux=   11.4,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h33m43.45', '-42.05.46.4'),
          rad=10727,
          flux=   13.5,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h33m43.89', '-45.27.18.4'),
          rad= 2877,
          flux=   83.5,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h33m44.44', '-45.48.34.1'),
          rad= 3752,
          flux=   39.9,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h33m54.40', '-48.20.52.1'),
          rad=12284,
          flux=   47.1,
          bmaj=  15.7,
          bmin= 0.0,
          bpa=  -2.3)

addsource(t, dir=dm.direction('J2000', '12h33m54.58', '-41.30.15.7'),
          rad=12834,
          flux=   98.0,
          bmaj=  56.6,
          bmin= 0.0,
          bpa= -21.3)

addsource(t, dir=dm.direction('J2000', '12h33m56.33', '-45.22.58.2'),
          rad= 2853,
          flux=   11.4,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h33m56.71', '-42.28.46.3'),
          rad= 9426,
          flux=   31.0,
          bmaj=  16.9,
          bmin= 0.0,
          bpa= -61.7)

addsource(t, dir=dm.direction('J2000', '12h34m02.35', '-44.04.44.6'),
          rad= 4207,
          flux=  125.7,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h34m03.26', '-46.53.02.3'),
          rad= 7240,
          flux=   61.0,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h34m07.75', '-43.31.14.9'),
          rad= 5952,
          flux=   10.6,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h34m08.00', '-48.00.18.3'),
          rad=11111,
          flux=   18.8,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h34m09.90', '-41.11.33.1'),
          rad=13966,
          flux=   32.8,
          bmaj=  16.7,
          bmin= 0.0,
          bpa=  45.5)

addsource(t, dir=dm.direction('J2000', '12h34m13.37', '-47.53.51.6'),
          rad=10750,
          flux=  355.9,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h34m16.46', '-42.48.37.3'),
          rad= 8353,
          flux=  164.9,
          bmaj=  65.4,
          bmin= 0.0,
          bpa= -15.2)

addsource(t, dir=dm.direction('J2000', '12h34m18.62', '-42.57.12.0'),
          rad= 7877,
          flux=  310.8,
          bmaj=  19.0,
          bmin= 0.0,
          bpa= -81.3)

addsource(t, dir=dm.direction('J2000', '12h34m20.86', '-49.14.27.4'),
          rad=15483,
          flux=  312.3,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h34m21.53', '-41.26.59.9'),
          rad=13087,
          flux=   35.6,
          bmaj=  21.8,
          bmin= 0.0,
          bpa=   1.2)

addsource(t, dir=dm.direction('J2000', '12h34m22.68', '-41.11.03.7'),
          rad=14023,
          flux=   18.7,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h34m28.16', '-41.09.54.1'),
          rad=14104,
          flux=  417.1,
          bmaj=  27.6,
          bmin= 0.0,
          bpa= -55.6)

addsource(t, dir=dm.direction('J2000', '12h34m32.42', '-41.09.14.9'),
          rad=14152,
          flux=  300.7,
          bmaj=  29.9,
          bmin= 0.0,
          bpa= -39.3)

addsource(t, dir=dm.direction('J2000', '12h34m32.45', '-45.52.03.2'),
          rad= 4240,
          flux=  108.1,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h34m33.15', '-47.49.10.6'),
          rad=10531,
          flux=   16.5,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h34m38.86', '-43.01.09.3'),
          rad= 7737,
          flux=   12.9,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h34m40.82', '-43.29.01.9'),
          rad= 6236,
          flux=   15.2,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h34m46.18', '-44.05.26.3'),
          rad= 4480,
          flux=   49.1,
          bmaj=  16.3,
          bmin= 0.0,
          bpa= -64.5)

addsource(t, dir=dm.direction('J2000', '12h34m46.84', '-47.45.32.8'),
          rad=10362,
          flux=   17.4,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h34m47.64', '-40.50.09.0'),
          rad=15306,
          flux=   40.2,
          bmaj=  63.2,
          bmin= 0.0,
          bpa= -57.5)

addsource(t, dir=dm.direction('J2000', '12h34m50.13', '-40.53.53.4'),
          rad=15092,
          flux=   11.1,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h34m51.24', '-45.33.43.9'),
          rad= 3680,
          flux=   70.3,
          bmaj=  52.4,
          bmin= 0.0,
          bpa=  13.1)

addsource(t, dir=dm.direction('J2000', '12h34m51.34', '-42.05.15.2'),
          rad=10947,
          flux=   58.8,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h34m52.08', '-46.30.10.5'),
          rad= 6213,
          flux=   19.4,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h34m53.12', '-44.26.58.1'),
          rad= 3699,
          flux=  312.5,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h34m54.89', '-40.34.51.4'),
          rad=16219,
          flux=   28.6,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h34m56.60', '-40.12.52.5'),
          rad=17514,
          flux=   44.4,
          bmaj=  21.0,
          bmin= 0.0,
          bpa=  55.7)

addsource(t, dir=dm.direction('J2000', '12h34m57.73', '-47.12.36.9'),
          rad= 8535,
          flux=   23.1,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h34m58.83', '-46.55.40.3'),
          rad= 7606,
          flux=   10.8,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h34m59.96', '-49.26.39.9'),
          rad=16272,
          flux=   17.2,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h35m06.43', '-40.15.10.5'),
          rad=17400,
          flux=   19.5,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h35m06.77', '-42.41.21.0'),
          rad= 8953,
          flux=   40.4,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h35m19.62', '-42.31.18.9'),
          rad= 9566,
          flux=   13.2,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h35m20.65', '-48.40.38.1'),
          rad=13630,
          flux=   17.3,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h35m20.69', '-40.59.09.9'),
          rad=14859,
          flux=   10.9,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h35m25.21', '-41.18.42.9'),
          rad=13735,
          flux=   27.9,
          bmaj=  33.0,
          bmin= 0.0,
          bpa= -15.8)

addsource(t, dir=dm.direction('J2000', '12h35m28.91', '-40.07.18.6'),
          rad=17910,
          flux=   60.9,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h35m30.12', '-46.41.31.8'),
          rad= 6999,
          flux=   25.3,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h35m34.34', '-43.50.01.9'),
          rad= 5518,
          flux=   12.5,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h35m35.50', '-42.06.18.3'),
          rad=11036,
          flux=   12.4,
          bmaj=  27.2,
          bmin= 0.0,
          bpa= -48.8)

addsource(t, dir=dm.direction('J2000', '12h35m35.72', '-41.37.07.6'),
          rad=12703,
          flux=  690.5,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h35m42.15', '-49.16.24.5'),
          rad=15759,
          flux=   38.9,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h35m42.96', '-41.18.52.7'),
          rad=13777,
          flux=   20.0,
          bmaj=  73.4,
          bmin= 0.0,
          bpa= -19.6)

addsource(t, dir=dm.direction('J2000', '12h35m46.92', '-42.19.00.2'),
          rad=10363,
          flux=   12.0,
          bmaj=  25.5,
          bmin= 0.0,
          bpa=   9.4)

addsource(t, dir=dm.direction('J2000', '12h35m48.47', '-41.57.01.5'),
          rad=11609,
          flux=   13.9,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h35m50.18', '-44.48.38.2'),
          rad= 3782,
          flux=   26.3,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h35m51.84', '-41.48.50.5'),
          rad=12086,
          flux=   11.3,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h35m54.98', '-42.55.23.8'),
          rad= 8399,
          flux=  112.2,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h35m55.33', '-43.54.47.6'),
          rad= 5456,
          flux=   15.6,
          bmaj=  25.1,
          bmin= 0.0,
          bpa=  77.4)

addsource(t, dir=dm.direction('J2000', '12h35m56.29', '-40.28.03.2'),
          rad=16763,
          flux=   44.0,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h35m57.39', '-46.04.54.9'),
          rad= 5409,
          flux=   27.1,
          bmaj=  23.5,
          bmin= 0.0,
          bpa=  77.6)

addsource(t, dir=dm.direction('J2000', '12h36m11.94', '-42.20.49.2'),
          rad=10363,
          flux=   15.8,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h36m12.09', '-46.49.37.2'),
          rad= 7636,
          flux=   87.6,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h36m20.66', '-43.19.45.8'),
          rad= 7275,
          flux=   17.0,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h36m23.93', '-46.02.58.7'),
          rad= 5527,
          flux=   28.2,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h36m26.18', '-44.55.42.9'),
          rad= 4106,
          flux=   13.8,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h36m29.59', '-47.23.33.2'),
          rad= 9512,
          flux=   18.2,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h36m30.91', '-47.01.02.3'),
          rad= 8324,
          flux=   19.3,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h36m32.38', '-47.16.00.9'),
          rad= 9120,
          flux=   14.9,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h36m34.00', '-41.39.43.1'),
          rad=12754,
          flux=   14.7,
          bmaj=  22.7,
          bmin= 0.0,
          bpa= -17.8)

addsource(t, dir=dm.direction('J2000', '12h36m34.10', '-49.02.01.2'),
          rad=15056,
          flux=   28.4,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h36m37.45', '-49.11.43.6'),
          rad=15623,
          flux=   57.4,
          bmaj=  78.7,
          bmin=  32.0,
          bpa=  21.3)

addsource(t, dir=dm.direction('J2000', '12h36m39.49', '-49.46.33.1'),
          rad=17642,
          flux=   22.8,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h36m40.40', '-45.03.45.1'),
          rad= 4250,
          flux=   20.4,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h36m43.34', '-46.54.17.5'),
          rad= 8042,
          flux=   35.0,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h36m43.51', '-49.31.21.9'),
          rad=16772,
          flux=   25.3,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h36m45.37', '-49.11.41.7'),
          rad=15643,
          flux=   43.3,
          bmaj=  52.8,
          bmin=  43.8,
          bpa=  32.0)

addsource(t, dir=dm.direction('J2000', '12h36m46.06', '-42.11.21.8'),
          rad=11032,
          flux=   20.0,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h36m46.20', '-46.04.43.3'),
          rad= 5769,
          flux=   26.8,
          bmaj=  32.5,
          bmin= 0.0,
          bpa= -27.7)

addsource(t, dir=dm.direction('J2000', '12h36m50.71', '-46.03.33.1'),
          rad= 5758,
          flux=   22.6,
          bmaj=  23.7,
          bmin= 0.0,
          bpa= -16.9)

addsource(t, dir=dm.direction('J2000', '12h36m50.87', '-45.25.37.5'),
          rad= 4605,
          flux=   23.6,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h36m51.02', '-47.00.28.8'),
          rad= 8399,
          flux=   27.8,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h36m51.67', '-47.05.48.8'),
          rad= 8677,
          flux=   19.7,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h36m53.25', '-41.57.50.2'),
          rad=11812,
          flux=   17.2,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h36m58.62', '-40.36.46.1'),
          rad=16433,
          flux=   11.7,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h36m59.54', '-42.36.42.0'),
          rad= 9720,
          flux=  194.4,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h37m00.41', '-43.30.00.6'),
          rad= 7038,
          flux=   19.2,
          bmaj=  20.2,
          bmin= 0.0,
          bpa= -44.5)

addsource(t, dir=dm.direction('J2000', '12h37m01.08', '-45.31.15.5'),
          rad= 4825,
          flux=  187.2,
          bmaj=  22.1,
          bmin= 0.0,
          bpa= -75.6)

addsource(t, dir=dm.direction('J2000', '12h37m02.07', '-42.59.23.1'),
          rad= 8548,
          flux=   13.0,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h37m02.95', '-49.13.43.8'),
          rad=15807,
          flux=   12.2,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h37m04.33', '-47.58.25.6'),
          rad=11561,
          flux=   12.3,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h37m04.56', '-45.17.27.6'),
          rad= 4612,
          flux=   20.1,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h37m07.58', '-49.13.02.9'),
          rad=15781,
          flux=   11.5,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h37m08.32', '-46.33.52.9'),
          rad= 7196,
          flux=   21.6,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h37m12.11', '-43.58.29.3'),
          rad= 5915,
          flux=   25.8,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h37m12.59', '-47.22.42.3'),
          rad= 9665,
          flux=   25.2,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h37m12.99', '-44.09.00.1'),
          rad= 5546,
          flux=   13.7,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h37m17.46', '-47.08.07.5'),
          rad= 8931,
          flux=   90.9,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h37m21.71', '-49.44.12.3'),
          rad=17610,
          flux=   10.3,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h37m23.44', '-43.23.32.7'),
          rad= 7497,
          flux=   31.5,
          bmaj=  32.1,
          bmin= 0.0,
          bpa= -42.9)

addsource(t, dir=dm.direction('J2000', '12h37m26.90', '-45.54.54.9'),
          rad= 5741,
          flux=   12.6,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h37m29.84', '-47.54.12.6'),
          rad=11433,
          flux=   13.6,
          bmaj=  26.7,
          bmin= 0.0,
          bpa= -10.5)

addsource(t, dir=dm.direction('J2000', '12h37m34.34', '-49.49.04.8'),
          rad=17923,
          flux=   13.5,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h37m34.52', '-41.35.29.4'),
          rad=13226,
          flux=   23.2,
          bmaj=  37.1,
          bmin= 0.0,
          bpa= -84.1)

addsource(t, dir=dm.direction('J2000', '12h37m39.14', '-45.16.03.5'),
          rad= 4953,
          flux=   14.0,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h37m42.77', '-47.39.04.2'),
          rad=10675,
          flux=   16.9,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h37m44.60', '-48.44.49.0'),
          rad=14293,
          flux=   12.8,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h37m45.06', '-47.36.59.0'),
          rad=10575,
          flux=   24.5,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h37m47.87', '-47.09.48.4'),
          rad= 9181,
          flux=   16.4,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h37m49.30', '-43.19.22.4'),
          rad= 7869,
          flux=   30.9,
          bmaj=  17.0,
          bmin= 0.0,
          bpa= -77.7)

addsource(t, dir=dm.direction('J2000', '12h37m50.86', '-46.31.27.1'),
          rad= 7373,
          flux=   23.8,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h37m53.19', '-41.44.34.3'),
          rad=12801,
          flux=   26.3,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h37m53.78', '-40.42.18.6'),
          rad=16297,
          flux=   14.6,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h37m57.52', '-46.05.29.4'),
          rad= 6371,
          flux=   22.5,
          bmaj=  27.2,
          bmin= 0.0,
          bpa=  -5.0)

addsource(t, dir=dm.direction('J2000', '12h37m57.92', '-45.45.38.2'),
          rad= 5731,
          flux=   13.4,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h38m00.15', '-40.24.20.0'),
          rad=17344,
          flux=   17.7,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h38m01.98', '-48.15.44.8'),
          rad=12742,
          flux=   13.7,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h38m08.91', '-46.45.60.0'),
          rad= 8153,
          flux=   14.4,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h38m11.69', '-41.53.04.0'),
          rad=12420,
          flux=   25.2,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h38m13.19', '-43.30.46.1'),
          rad= 7530,
          flux=   13.9,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h38m18.21', '-48.04.14.2'),
          rad=12183,
          flux=   14.9,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h38m18.67', '-49.40.29.0'),
          rad=17553,
          flux=   14.9,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h38m19.05', '-41.24.01.4'),
          rad=14048,
          flux=   43.7,
          bmaj=  49.7,
          bmin= 0.0,
          bpa=  21.6)

addsource(t, dir=dm.direction('J2000', '12h38m19.10', '-47.17.18.5'),
          rad= 9731,
          flux=   35.8,
          bmaj=  22.8,
          bmin= 0.0,
          bpa=  86.8)

addsource(t, dir=dm.direction('J2000', '12h38m23.77', '-49.06.60.0'),
          rad=15672,
          flux=   21.6,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h38m24.96', '-42.28.37.0'),
          rad=10599,
          flux=   77.3,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h38m26.27', '-44.31.38.1'),
          rad= 5653,
          flux=   12.4,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h38m29.45', '-45.46.54.0'),
          rad= 6058,
          flux=   72.2,
          bmaj=  39.3,
          bmin= 0.0,
          bpa= -33.6)

addsource(t, dir=dm.direction('J2000', '12h38m30.50', '-40.44.34.2'),
          rad=16302,
          flux=   23.0,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h38m32.36', '-47.57.48.6'),
          rad=11901,
          flux=   28.9,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h38m33.81', '-44.13.46.7'),
          rad= 6146,
          flux=  239.4,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h38m36.59', '-44.47.14.4'),
          rad= 5542,
          flux=   10.7,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h38m37.67', '-45.39.15.9'),
          rad= 5945,
          flux=   82.0,
          bmaj=  15.0,
          bmin= 0.0,
          bpa= -49.7)

addsource(t, dir=dm.direction('J2000', '12h38m38.39', '-44.10.56.8'),
          rad= 6270,
          flux=   82.8,
          bmaj=  80.6,
          bmin=  35.8,
          bpa= -59.7)

addsource(t, dir=dm.direction('J2000', '12h38m39.73', '-44.42.13.6'),
          rad= 5628,
          flux=   11.6,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h38m43.68', '-45.54.52.8'),
          rad= 6417,
          flux=   71.3,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h38m47.82', '-44.09.35.3'),
          rad= 6398,
          flux=   23.0,
          bmaj=  59.6,
          bmin=  43.9,
          bpa= -42.4)

addsource(t, dir=dm.direction('J2000', '12h38m48.30', '-47.26.56.6'),
          rad=10377,
          flux=   22.4,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h38m49.00', '-46.07.52.2'),
          rad= 6886,
          flux=   34.7,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h38m49.51', '-42.30.26.5'),
          rad=10645,
          flux=   10.6,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h38m49.75', '-45.11.04.6'),
          rad= 5648,
          flux=   15.8,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h38m55.34', '-44.08.40.9'),
          rad= 6495,
          flux=  142.2,
          bmaj=  69.5,
          bmin=  24.4,
          bpa= -56.6)

addsource(t, dir=dm.direction('J2000', '12h38m57.84', '-48.38.25.1'),
          rad=14208,
          flux=   23.1,
          bmaj=  61.7,
          bmin= 0.0,
          bpa=  58.7)

addsource(t, dir=dm.direction('J2000', '12h39m02.35', '-42.52.05.5'),
          rad= 9651,
          flux=   29.1,
          bmaj=  61.7,
          bmin= 0.0,
          bpa=  75.2)

addsource(t, dir=dm.direction('J2000', '12h39m03.75', '-44.47.03.5'),
          rad= 5829,
          flux=   12.2,
          bmaj=  29.1,
          bmin= 0.0,
          bpa=  83.4)

addsource(t, dir=dm.direction('J2000', '12h39m03.92', '-41.49.25.0'),
          rad=12870,
          flux=   38.8,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h39m03.92', '-44.27.01.8'),
          rad= 6124,
          flux=   19.1,
          bmaj=  32.5,
          bmin= 0.0,
          bpa=  67.1)

addsource(t, dir=dm.direction('J2000', '12h39m04.91', '-40.57.46.3'),
          rad=15699,
          flux=   13.5,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h39m05.82', '-41.00.14.5'),
          rad=15565,
          flux=   18.3,
          bmaj=  24.0,
          bmin= 0.0,
          bpa= -12.2)

addsource(t, dir=dm.direction('J2000', '12h39m09.50', '-44.23.35.7'),
          rad= 6252,
          flux=   43.0,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h39m09.96', '-44.00.29.3'),
          rad= 6881,
          flux=   18.9,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h39m10.37', '-48.33.51.4'),
          rad=14009,
          flux=   80.6,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h39m10.88', '-44.37.41.6'),
          rad= 6012,
          flux=   11.0,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h39m14.53', '-41.28.37.0'),
          rad=14044,
          flux=   19.2,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h39m15.03', '-40.49.56.6'),
          rad=16176,
          flux=   69.7,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h39m15.89', '-41.37.02.6'),
          rad=13594,
          flux=  140.1,
          bmaj=  20.9,
          bmin= 0.0,
          bpa= -69.3)

addsource(t, dir=dm.direction('J2000', '12h39m19.27', '-40.26.19.8'),
          rad=17516,
          flux=   28.2,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h39m25.52', '-47.28.17.1'),
          rad=10652,
          flux=  122.6,
          bmaj=  48.4,
          bmin= 0.0,
          bpa= -29.5)

addsource(t, dir=dm.direction('J2000', '12h39m26.33', '-43.52.18.3'),
          rad= 7298,
          flux=   44.3,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h39m27.33', '-43.50.07.9'),
          rad= 7382,
          flux=   52.9,
          bmaj=  15.2,
          bmin= 0.0,
          bpa=   1.5)

addsource(t, dir=dm.direction('J2000', '12h39m28.32', '-48.26.49.2'),
          rad=13705,
          flux=   11.5,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h39m37.97', '-40.32.50.9'),
          rad=17224,
          flux=   19.3,
          bmaj=  36.0,
          bmin= 0.0,
          bpa=  37.6)

addsource(t, dir=dm.direction('J2000', '12h39m40.30', '-41.04.12.5'),
          rad=15496,
          flux=   11.8,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h39m42.03', '-44.15.54.8'),
          rad= 6751,
          flux=   72.9,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h39m43.06', '-44.36.52.2'),
          rad= 6357,
          flux=   20.4,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h39m46.84', '-41.04.26.4'),
          rad=15513,
          flux=   14.7,
          bmaj=  44.6,
          bmin= 0.0,
          bpa= -73.5)

addsource(t, dir=dm.direction('J2000', '12h39m47.63', '-42.14.53.2'),
          rad=11776,
          flux=   20.5,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h39m52.97', '-42.26.57.5'),
          rad=11202,
          flux=   33.0,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h39m55.70', '-43.16.08.9'),
          rad= 8938,
          flux=   17.7,
          bmaj=  38.1,
          bmin=  29.9,
          bpa= -67.4)

addsource(t, dir=dm.direction('J2000', '12h39m58.14', '-41.24.26.2'),
          rad=14480,
          flux=   24.6,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h40m00.21', '-46.27.24.4'),
          rad= 8183,
          flux=   64.2,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h40m01.48', '-46.34.33.5'),
          rad= 8468,
          flux=   16.0,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h40m01.56', '-44.46.16.0'),
          rad= 6445,
          flux=   47.4,
          bmaj=  36.5,
          bmin= 0.0,
          bpa=  58.3)

addsource(t, dir=dm.direction('J2000', '12h40m03.07', '-41.59.04.1'),
          rad=12676,
          flux=   96.8,
          bmaj=  20.4,
          bmin= 0.0,
          bpa=  14.3)

addsource(t, dir=dm.direction('J2000', '12h40m05.27', '-41.14.41.8'),
          rad=15040,
          flux=   23.8,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h40m09.24', '-47.39.44.5'),
          rad=11468,
          flux=   22.9,
          bmaj=  60.6,
          bmin= 0.0,
          bpa= -47.6)

addsource(t, dir=dm.direction('J2000', '12h40m11.44', '-48.47.11.1'),
          rad=14988,
          flux=   13.4,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h40m22.07', '-46.56.37.7'),
          rad= 9536,
          flux=  142.8,
          bmaj=  42.1,
          bmin=  15.3,
          bpa= -21.6)

addsource(t, dir=dm.direction('J2000', '12h40m22.35', '-44.50.46.3'),
          rad= 6632,
          flux=   13.4,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h40m23.05', '-40.58.50.8'),
          rad=15985,
          flux=   82.6,
          bmaj=  15.8,
          bmin= 0.0,
          bpa= -53.5)

addsource(t, dir=dm.direction('J2000', '12h40m25.05', '-47.01.59.8'),
          rad= 9792,
          flux=   37.2,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h40m30.73', '-46.04.40.6'),
          rad= 7677,
          flux=   36.3,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h40m31.54', '-48.59.04.3'),
          rad=15715,
          flux=   84.3,
          bmaj=  23.2,
          bmin= 0.0,
          bpa=  14.4)

addsource(t, dir=dm.direction('J2000', '12h40m32.50', '-43.39.11.4'),
          rad= 8338,
          flux=   49.8,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h40m32.54', '-49.00.26.2'),
          rad=15793,
          flux=   37.5,
          bmaj=  28.0,
          bmin= 0.0,
          bpa=  18.3)

addsource(t, dir=dm.direction('J2000', '12h40m33.58', '-42.46.51.2'),
          rad=10517,
          flux=   33.8,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h40m35.51', '-41.13.22.3'),
          rad=15259,
          flux=   27.5,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h40m40.70', '-45.51.55.5'),
          rad= 7427,
          flux=   13.4,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h40m41.72', '-46.27.34.2'),
          rad= 8526,
          flux=   98.8,
          bmaj=  36.2,
          bmin= 0.0,
          bpa=  35.2)

addsource(t, dir=dm.direction('J2000', '12h40m42.50', '-42.16.55.4'),
          rad=12008,
          flux=   23.1,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h40m46.68', '-44.47.50.8'),
          rad= 6908,
          flux=   15.4,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h40m49.18', '-42.03.37.5'),
          rad=12711,
          flux=   15.0,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h40m51.24', '-49.08.03.4'),
          rad=16283,
          flux=   65.1,
          bmaj=  21.5,
          bmin= 0.0,
          bpa=  71.1)

addsource(t, dir=dm.direction('J2000', '12h40m53.87', '-44.10.29.7'),
          rad= 7588,
          flux=  121.2,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h40m59.81', '-40.48.45.7'),
          rad=16706,
          flux=   10.9,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h40m59.89', '-44.07.31.9'),
          rad= 7721,
          flux=  130.0,
          bmaj=  65.0,
          bmin= 0.0,
          bpa=  70.7)

addsource(t, dir=dm.direction('J2000', '12h41m00.67', '-46.32.34.9'),
          rad= 8864,
          flux=   17.6,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h41m03.05', '-47.45.38.9'),
          rad=12068,
          flux=   13.7,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h41m03.21', '-41.30.42.9'),
          rad=14483,
          flux=   67.0,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h41m03.84', '-45.53.06.2'),
          rad= 7676,
          flux=   16.8,
          bmaj=  25.9,
          bmin= 0.0,
          bpa=   4.7)

addsource(t, dir=dm.direction('J2000', '12h41m04.44', '-47.36.33.3'),
          rad=11639,
          flux=   68.2,
          bmaj=  34.0,
          bmin= 0.0,
          bpa=  38.5)

addsource(t, dir=dm.direction('J2000', '12h41m06.62', '-48.27.16.6'),
          rad=14187,
          flux=   87.8,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h41m09.08', '-49.32.33.3'),
          rad=17690,
          flux=  111.5,
          bmaj=  32.8,
          bmin= 0.0,
          bpa= -30.6)

addsource(t, dir=dm.direction('J2000', '12h41m09.55', '-47.06.36.5'),
          rad=10304,
          flux=   22.9,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h41m09.84', '-46.34.45.9'),
          rad= 9019,
          flux=  270.4,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h41m11.45', '-49.00.04.3'),
          rad=15939,
          flux=   20.0,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h41m12.42', '-41.26.41.5'),
          rad=14744,
          flux=   13.4,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h41m12.45', '-40.45.19.1'),
          rad=16953,
          flux=   11.7,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h41m13.38', '-46.50.24.0'),
          rad= 9652,
          flux=   67.1,
          bmaj=  23.3,
          bmin= 0.0,
          bpa= -51.6)

addsource(t, dir=dm.direction('J2000', '12h41m13.58', '-44.24.56.1'),
          rad= 7481,
          flux=   27.4,
          bmaj=  20.2,
          bmin= 0.0,
          bpa= -46.0)

addsource(t, dir=dm.direction('J2000', '12h41m14.94', '-43.56.26.4'),
          rad= 8167,
          flux=   15.4,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h41m21.75', '-42.49.10.2'),
          rad=10759,
          flux=   15.1,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h41m26.16', '-43.51.42.9'),
          rad= 8411,
          flux=   27.8,
          bmaj=  24.6,
          bmin= 0.0,
          bpa= -53.1)

addsource(t, dir=dm.direction('J2000', '12h41m27.35', '-46.00.52.3'),
          rad= 8093,
          flux=   23.0,
          bmaj=  31.6,
          bmin= 0.0,
          bpa= -72.8)

addsource(t, dir=dm.direction('J2000', '12h41m29.93', '-42.39.39.1'),
          rad=11247,
          flux=   98.2,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h41m31.38', '-45.23.24.2'),
          rad= 7440,
          flux=   96.3,
          bmaj=  72.6,
          bmin= 0.0,
          bpa= -49.5)

addsource(t, dir=dm.direction('J2000', '12h41m31.48', '-43.40.16.6'),
          rad= 8824,
          flux=  220.7,
          bmaj=  14.6,
          bmin= 0.0,
          bpa=   0.9)

addsource(t, dir=dm.direction('J2000', '12h41m34.57', '-43.58.49.7'),
          rad= 8286,
          flux=   67.6,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h41m35.33', '-40.55.39.2'),
          rad=16508,
          flux=   26.4,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h41m37.49', '-45.06.29.7'),
          rad= 7399,
          flux=   13.4,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h41m39.24', '-44.19.57.4'),
          rad= 7835,
          flux=   10.7,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h41m40.60', '-44.08.37.5'),
          rad= 8093,
          flux=   15.5,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h41m40.74', '-49.34.12.0'),
          rad=17904,
          flux=  140.2,
          bmaj=  65.1,
          bmin=  15.8,
          bpa=  62.7)

addsource(t, dir=dm.direction('J2000', '12h41m47.19', '-44.54.31.4'),
          rad= 7512,
          flux=   16.6,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h41m47.69', '-43.02.33.8'),
          rad=10383,
          flux=   26.3,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h41m55.30', '-48.56.21.1'),
          rad=15941,
          flux=   26.3,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h41m58.44', '-45.00.54.4'),
          rad= 7617,
          flux=   15.4,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h41m58.65', '-43.08.33.0'),
          rad=10227,
          flux=   31.9,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h41m58.93', '-47.22.59.6'),
          rad=11365,
          flux=   17.2,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h41m59.06', '-43.16.25.3'),
          rad= 9922,
          flux=   49.1,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h41m59.45', '-42.36.52.4'),
          rad=11585,
          flux=   85.7,
          bmaj=  15.5,
          bmin= 0.0,
          bpa= -86.3)

addsource(t, dir=dm.direction('J2000', '12h42m01.70', '-41.03.37.8'),
          rad=16221,
          flux=   16.8,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h42m03.69', '-43.27.59.0'),
          rad= 9534,
          flux=   75.2,
          bmaj=  26.0,
          bmin= 0.0,
          bpa=  51.6)

addsource(t, dir=dm.direction('J2000', '12h42m06.50', '-49.31.15.2'),
          rad=17851,
          flux=   13.6,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h42m11.97', '-42.16.25.3'),
          rad=12618,
          flux=   50.9,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h42m15.66', '-44.07.15.5'),
          rad= 8472,
          flux=   38.1,
          bmaj=  40.8,
          bmin=  16.9,
          bpa= -55.8)

addsource(t, dir=dm.direction('J2000', '12h42m16.67', '-45.23.54.2'),
          rad= 7915,
          flux=   33.4,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h42m18.79', '-43.55.59.9'),
          rad= 8788,
          flux=  485.5,
          bmaj=  37.4,
          bmin= 0.0,
          bpa=  36.2)

addsource(t, dir=dm.direction('J2000', '12h42m21.74', '-43.56.47.4'),
          rad= 8795,
          flux=   23.1,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h42m22.87', '-44.59.59.2'),
          rad= 7877,
          flux=   31.6,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h42m26.92', '-48.39.48.2'),
          rad=15237,
          flux=   10.6,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h42m28.58', '-43.02.07.2'),
          rad=10728,
          flux=   25.2,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h42m30.89', '-44.02.16.3'),
          rad= 8743,
          flux=  257.0,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h42m31.09', '-45.02.53.4'),
          rad= 7963,
          flux=   16.3,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h42m31.52', '-42.40.10.7'),
          rad=11675,
          flux=   37.1,
          bmaj=  16.3,
          bmin= 0.0,
          bpa=   8.9)

addsource(t, dir=dm.direction('J2000', '12h42m34.21', '-43.02.20.6'),
          rad=10765,
          flux=   27.4,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h42m37.75', '-42.44.05.4'),
          rad=11553,
          flux=  273.2,
          bmaj=  30.7,
          bmin= 0.0,
          bpa= -77.9)

addsource(t, dir=dm.direction('J2000', '12h42m38.69', '-46.23.15.5'),
          rad= 9384,
          flux=   27.8,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h42m39.13', '-48.43.38.3'),
          rad=15496,
          flux=   66.4,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h42m41.79', '-42.34.45.0'),
          rad=11991,
          flux=  116.7,
          bmaj=  44.2,
          bmin= 0.0,
          bpa= -18.1)

addsource(t, dir=dm.direction('J2000', '12h42m44.92', '-43.05.53.7'),
          rad=10712,
          flux=   53.7,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h42m47.79', '-40.41.47.3'),
          rad=17619,
          flux=   43.3,
          bmaj=  21.5,
          bmin= 0.0,
          bpa=  59.1)

addsource(t, dir=dm.direction('J2000', '12h42m47.90', '-45.02.21.9'),
          rad= 8141,
          flux=   24.5,
          bmaj=  25.5,
          bmin= 0.0,
          bpa=  39.5)

addsource(t, dir=dm.direction('J2000', '12h42m48.19', '-42.42.06.9'),
          rad=11718,
          flux=   29.9,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h42m48.47', '-40.56.20.8'),
          rad=16856,
          flux=   13.8,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h42m48.57', '-49.23.22.1'),
          rad=17612,
          flux=   15.4,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h42m49.49', '-44.22.01.1'),
          rad= 8514,
          flux=   15.3,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h42m53.32', '-45.27.43.1'),
          rad= 8334,
          flux=  152.5,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h42m53.73', '-43.29.02.0'),
          rad= 9941,
          flux=   27.6,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h42m55.96', '-44.36.21.5'),
          rad= 8377,
          flux=   10.9,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h42m56.54', '-44.51.46.5'),
          rad= 8258,
          flux=   26.6,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h42m57.61', '-45.26.28.4'),
          rad= 8365,
          flux=   13.7,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h43m00.30', '-49.05.43.9'),
          rad=16740,
          flux=   40.9,
          bmaj=  34.8,
          bmin= 0.0,
          bpa= -51.9)

addsource(t, dir=dm.direction('J2000', '12h43m02.46', '-45.05.37.0'),
          rad= 8296,
          flux=   21.3,
          bmaj=  21.3,
          bmin= 0.0,
          bpa=  -8.1)

addsource(t, dir=dm.direction('J2000', '12h43m07.76', '-41.37.42.0'),
          rad=14859,
          flux=   11.6,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h43m08.16', '-43.12.40.6'),
          rad=10649,
          flux=   61.9,
          bmaj=  27.6,
          bmin= 0.0,
          bpa=  67.3)

addsource(t, dir=dm.direction('J2000', '12h43m14.88', '-44.49.02.7'),
          rad= 8467,
          flux=   20.5,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h43m15.12', '-43.13.01.0'),
          rad=10696,
          flux=  192.2,
          bmaj=  15.9,
          bmin= 0.0,
          bpa=  74.3)

addsource(t, dir=dm.direction('J2000', '12h43m18.02', '-46.02.17.7'),
          rad= 9178,
          flux=   12.6,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h43m19.82', '-44.00.26.4'),
          rad= 9269,
          flux=   52.7,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h43m21.45', '-48.17.44.6'),
          rad=14438,
          flux=   95.1,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h43m24.87', '-41.45.30.4'),
          rad=14586,
          flux=   16.4,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h43m31.54', '-41.54.10.4'),
          rad=14213,
          flux=   80.1,
          bmaj=  16.1,
          bmin= 0.0,
          bpa=   1.3)

addsource(t, dir=dm.direction('J2000', '12h43m34.46', '-41.28.15.1'),
          rad=15494,
          flux=   14.0,
          bmaj=  37.2,
          bmin= 0.0,
          bpa=   9.3)

addsource(t, dir=dm.direction('J2000', '12h43m36.86', '-41.29.10.6'),
          rad=15464,
          flux=   18.1,
          bmaj=  23.7,
          bmin= 0.0,
          bpa=  18.2)

addsource(t, dir=dm.direction('J2000', '12h43m39.89', '-47.55.00.2'),
          rad=13479,
          flux=   16.9,
          bmaj=  22.1,
          bmin= 0.0,
          bpa=  62.9)

addsource(t, dir=dm.direction('J2000', '12h43m49.85', '-43.38.16.4'),
          rad=10161,
          flux=   41.9,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h43m55.75', '-40.43.59.1'),
          rad=17873,
          flux=  186.3,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h43m58.17', '-47.29.54.6'),
          rad=12499,
          flux=   91.8,
          bmaj=  17.6,
          bmin= 0.0,
          bpa= -12.6)

addsource(t, dir=dm.direction('J2000', '12h44m06.04', '-42.39.01.2'),
          rad=12455,
          flux=   12.7,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h44m16.27', '-42.15.16.1'),
          rad=13557,
          flux=  273.7,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h44m17.78', '-45.19.54.1'),
          rad= 9146,
          flux=   89.3,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h44m18.66', '-48.40.34.5'),
          rad=15879,
          flux=   11.7,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h44m19.96', '-40.51.36.8'),
          rad=17618,
          flux=  361.0,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h44m26.02', '-45.45.26.4'),
          rad= 9519,
          flux=  171.2,
          bmaj=  16.5,
          bmin= 0.0,
          bpa=  36.5)

addsource(t, dir=dm.direction('J2000', '12h44m29.86', '-47.59.55.1'),
          rad=14029,
          flux=   11.6,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h44m33.46', '-46.09.08.2'),
          rad=10060,
          flux=  251.9,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h44m34.40', '-41.35.53.9'),
          rad=15509,
          flux=   97.7,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h44m34.69', '-42.54.22.7'),
          rad=12075,
          flux=   53.4,
          bmaj=  19.4,
          bmin= 0.0,
          bpa=  65.9)

addsource(t, dir=dm.direction('J2000', '12h44m35.69', '-43.31.05.1'),
          rad=10808,
          flux=   46.3,
          bmaj=  28.8,
          bmin= 0.0,
          bpa=  59.7)

addsource(t, dir=dm.direction('J2000', '12h44m39.69', '-46.18.57.5'),
          rad=10363,
          flux=   44.9,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h44m40.42', '-42.22.22.0'),
          rad=13428,
          flux=   57.5,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h44m42.61', '-45.52.46.3'),
          rad= 9809,
          flux=   11.7,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h44m44.27', '-48.48.56.1'),
          rad=16434,
          flux=   59.1,
          bmaj=  28.1,
          bmin= 0.0,
          bpa= -54.4)

addsource(t, dir=dm.direction('J2000', '12h44m44.65', '-45.16.08.3'),
          rad= 9407,
          flux=   30.4,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h44m45.07', '-44.15.34.3'),
          rad= 9812,
          flux=   11.1,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h44m47.21', '-44.31.35.0'),
          rad= 9597,
          flux=  106.7,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h44m48.58', '-41.36.35.0'),
          rad=15571,
          flux=   22.2,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h44m49.18', '-40.48.06.5'),
          rad=17969,
          flux= 1341.8,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h44m55.58', '-42.57.06.6'),
          rad=12149,
          flux=   18.8,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h44m57.34', '-45.33.49.0'),
          rad= 9681,
          flux=   11.0,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h44m57.90', '-44.47.23.7'),
          rad= 9567,
          flux=   22.8,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h44m59.29', '-49.04.52.6'),
          rad=17306,
          flux=  210.6,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h45m00.51', '-43.36.00.2'),
          rad=10896,
          flux=  159.4,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h45m02.11', '-42.36.53.0'),
          rad=12993,
          flux=   47.5,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h45m04.17', '-47.53.11.5'),
          rad=13962,
          flux=   43.5,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h45m12.75', '-41.32.18.6'),
          rad=15939,
          flux=   34.1,
          bmaj=  55.7,
          bmin= 0.0,
          bpa=  61.9)

addsource(t, dir=dm.direction('J2000', '12h45m15.55', '-44.01.47.3'),
          rad=10391,
          flux=   11.4,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h45m15.69', '-43.40.20.3'),
          rad=10919,
          flux=   20.2,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h45m17.28', '-44.50.49.7'),
          rad= 9753,
          flux=   14.1,
          bmaj=  25.3,
          bmin= 0.0,
          bpa=  24.2)

addsource(t, dir=dm.direction('J2000', '12h45m17.81', '-43.38.56.4'),
          rad=10978,
          flux=   78.9,
          bmaj=  27.0,
          bmin= 0.0,
          bpa=  56.3)

addsource(t, dir=dm.direction('J2000', '12h45m23.29', '-45.25.07.6'),
          rad= 9868,
          flux=   10.9,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h45m26.80', '-43.24.56.2'),
          rad=11474,
          flux=   32.6,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h45m29.27', '-43.58.02.8'),
          rad=10611,
          flux=   22.4,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h45m29.56', '-43.21.35.2'),
          rad=11605,
          flux=   13.6,
          bmaj=  29.0,
          bmin= 0.0,
          bpa= -18.3)

addsource(t, dir=dm.direction('J2000', '12h45m30.17', '-44.46.52.2'),
          rad= 9911,
          flux=  185.8,
          bmaj=  54.1,
          bmin= 0.0,
          bpa=   3.8)

addsource(t, dir=dm.direction('J2000', '12h45m32.23', '-45.17.08.8'),
          rad= 9912,
          flux=   10.9,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h45m33.92', '-45.44.22.9'),
          rad=10190,
          flux=   33.3,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h45m34.26', '-41.39.57.5'),
          rad=15728,
          flux=   14.9,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h45m36.81', '-45.25.05.7'),
          rad=10009,
          flux=   45.9,
          bmaj=  35.6,
          bmin= 0.0,
          bpa=  25.9)

addsource(t, dir=dm.direction('J2000', '12h45m38.15', '-43.23.51.3'),
          rad=11614,
          flux=   11.2,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h45m43.56', '-48.59.24.7'),
          rad=17281,
          flux=  206.1,
          bmaj=  23.2,
          bmin= 0.0,
          bpa= -13.0)

addsource(t, dir=dm.direction('J2000', '12h45m52.92', '-41.50.28.9'),
          rad=15379,
          flux=   13.9,
          bmaj=  41.3,
          bmin= 0.0,
          bpa=   6.0)

addsource(t, dir=dm.direction('J2000', '12h45m58.58', '-41.04.30.2'),
          rad=17583,
          flux=   13.4,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h45m59.54', '-41.38.42.0'),
          rad=15965,
          flux=   12.3,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h45m59.94', '-43.46.27.0'),
          rad=11189,
          flux=   75.9,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h46m01.33', '-44.58.54.4'),
          rad=10193,
          flux=   11.0,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h46m06.04', '-48.01.15.7'),
          rad=14739,
          flux=   13.5,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h46m08.07', '-45.43.33.5'),
          rad=10526,
          flux=   17.7,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h46m09.14', '-48.45.17.7'),
          rad=16753,
          flux=   14.9,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h46m20.18', '-48.05.08.6'),
          rad=15005,
          flux=   42.8,
          bmaj=  67.2,
          bmin= 0.0,
          bpa=  62.5)

addsource(t, dir=dm.direction('J2000', '12h46m24.52', '-48.02.34.5'),
          rad=14924,
          flux=   16.4,
          bmaj=  21.3,
          bmin= 0.0,
          bpa=  40.5)

addsource(t, dir=dm.direction('J2000', '12h46m24.68', '-47.18.28.5'),
          rad=13168,
          flux=   50.2,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h46m29.47', '-44.41.28.6'),
          rad=10576,
          flux=   17.9,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h46m31.64', '-44.19.17.2'),
          rad=10852,
          flux=   19.1,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h46m38.92', '-43.18.06.9'),
          rad=12358,
          flux=   14.3,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h46m39.56', '-43.23.30.7'),
          rad=12200,
          flux=   12.3,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h46m41.10', '-47.40.34.7'),
          rad=14138,
          flux=   27.3,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h46m49.54', '-42.07.08.2'),
          rad=15082,
          flux=  204.1,
          bmaj=  17.7,
          bmin= 0.0,
          bpa=  31.7)

addsource(t, dir=dm.direction('J2000', '12h46m53.32', '-45.55.11.3'),
          rad=11157,
          flux=   23.2,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h46m56.70', '-42.16.00.4'),
          rad=14770,
          flux=   31.3,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h47m02.01', '-45.52.55.4'),
          rad=11208,
          flux=   19.7,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h47m02.47', '-42.18.43.7'),
          rad=14706,
          flux=   16.2,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h47m07.40', '-49.00.18.2'),
          rad=17813,
          flux= 1204.6,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h47m09.17', '-43.49.31.7'),
          rad=11802,
          flux=   27.0,
          bmaj=  71.0,
          bmin=  23.8,
          bpa= -41.5)

addsource(t, dir=dm.direction('J2000', '12h47m12.33', '-42.39.39.7'),
          rad=13975,
          flux=   15.6,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h47m13.89', '-43.28.32.0'),
          rad=12383,
          flux=   29.6,
          bmaj=  21.2,
          bmin= 0.0,
          bpa= -24.8)

addsource(t, dir=dm.direction('J2000', '12h47m18.16', '-43.50.33.4'),
          rad=11868,
          flux=   11.6,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h47m21.03', '-43.55.59.4'),
          rad=11779,
          flux=  106.7,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h47m21.94', '-49.01.56.3'),
          rad=17978,
          flux=   10.5,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h47m22.57', '-48.16.00.1'),
          rad=15904,
          flux=  143.7,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h47m22.74', '-47.32.22.4'),
          rad=14144,
          flux=   17.7,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h47m26.70', '-46.46.07.7'),
          rad=12638,
          flux=  134.6,
          bmaj=  19.7,
          bmin= 0.0,
          bpa= -48.2)

addsource(t, dir=dm.direction('J2000', '12h47m30.27', '-48.07.01.1'),
          rad=15577,
          flux=   22.1,
          bmaj=  74.4,
          bmin= 0.0,
          bpa=  21.8)

addsource(t, dir=dm.direction('J2000', '12h47m30.65', '-46.00.43.4'),
          rad=11622,
          flux=   21.4,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h47m30.96', '-41.50.14.9'),
          rad=16126,
          flux=   53.1,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h47m32.09', '-41.12.09.0'),
          rad=17851,
          flux=   11.4,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h47m34.48', '-41.37.27.4'),
          rad=16714,
          flux=   15.2,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h47m34.69', '-43.09.09.0'),
          rad=13157,
          flux=   12.6,
          bmaj=  27.0,
          bmin= 0.0,
          bpa=  53.0)

addsource(t, dir=dm.direction('J2000', '12h47m35.15', '-41.47.33.3'),
          rad=16275,
          flux=   20.0,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h47m43.83', '-48.51.58.0'),
          rad=17648,
          flux=   14.4,
          bmaj=  24.2,
          bmin= 0.0,
          bpa=  -1.1)

addsource(t, dir=dm.direction('J2000', '12h47m43.99', '-43.20.16.6'),
          rad=12906,
          flux=   31.0,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h47m48.50', '-45.39.17.8'),
          rad=11505,
          flux=  101.9,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h47m54.66', '-46.10.45.3'),
          rad=12044,
          flux=   53.8,
          bmaj=  24.6,
          bmin= 0.0,
          bpa= -57.8)

addsource(t, dir=dm.direction('J2000', '12h47m58.01', '-45.49.48.5'),
          rad=11730,
          flux=   25.3,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h47m58.47', '-42.51.10.3'),
          rad=13968,
          flux=   67.1,
          bmaj=  64.7,
          bmin= 0.0,
          bpa= -37.8)

addsource(t, dir=dm.direction('J2000', '12h47m59.11', '-46.10.16.1'),
          rad=12078,
          flux=   47.2,
          bmaj=  26.7,
          bmin= 0.0,
          bpa= -57.9)

addsource(t, dir=dm.direction('J2000', '12h48m06.82', '-47.16.59.3'),
          rad=13955,
          flux=   16.7,
          bmaj=  52.0,
          bmin= 0.0,
          bpa= -39.1)

addsource(t, dir=dm.direction('J2000', '12h48m14.79', '-47.23.39.4'),
          rad=14250,
          flux=   25.5,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h48m16.08', '-42.08.27.9'),
          rad=15723,
          flux=   13.9,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h48m23.18', '-41.54.56.6'),
          rad=16335,
          flux=   32.4,
          bmaj=  64.0,
          bmin= 0.0,
          bpa= -11.0)

addsource(t, dir=dm.direction('J2000', '12h48m26.52', '-46.42.05.9'),
          rad=13072,
          flux=  108.5,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h48m31.68', '-48.51.48.2'),
          rad=17945,
          flux=   15.5,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h48m34.54', '-47.11.52.5'),
          rad=14021,
          flux=   26.1,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h48m36.77', '-45.10.22.2'),
          rad=11835,
          flux=   13.3,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h48m47.64', '-47.59.55.4'),
          rad=15857,
          flux=   54.6,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h48m48.81', '-41.35.13.8'),
          rad=17376,
          flux=   19.8,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h48m50.95', '-48.39.53.4'),
          rad=17543,
          flux=  106.5,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h48m52.24', '-46.26.50.9'),
          rad=12940,
          flux=   19.2,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h48m57.24', '-47.03.44.5'),
          rad=13964,
          flux=  338.7,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h49m00.45', '-47.14.27.5'),
          rad=14326,
          flux=   91.4,
          bmaj=  14.6,
          bmin= 0.0,
          bpa= -10.7)

addsource(t, dir=dm.direction('J2000', '12h49m05.18', '-43.10.17.9'),
          rad=13971,
          flux=   12.7,
          bmaj=  28.3,
          bmin= 0.0,
          bpa= -54.3)

addsource(t, dir=dm.direction('J2000', '12h49m10.28', '-41.39.04.2'),
          rad=17377,
          flux=   25.7,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h49m11.95', '-47.48.36.8'),
          rad=15610,
          flux=   22.6,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h49m12.15', '-44.17.52.1'),
          rad=12543,
          flux=  467.3,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h49m13.52', '-45.45.30.8'),
          rad=12447,
          flux=   59.1,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h49m16.60', '-42.01.49.3'),
          rad=16487,
          flux=  179.6,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h49m18.98', '-43.52.59.5'),
          rad=13036,
          flux=   12.0,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h49m23.00', '-44.44.45.1'),
          rad=12387,
          flux=  681.2,
          bmaj=  24.1,
          bmin= 0.0,
          bpa=  47.4)

addsource(t, dir=dm.direction('J2000', '12h49m23.37', '-46.07.07.4'),
          rad=12854,
          flux=   29.7,
          bmaj=  54.8,
          bmin= 0.0,
          bpa=   0.7)

addsource(t, dir=dm.direction('J2000', '12h49m25.69', '-48.34.07.4'),
          rad=17533,
          flux=   19.1,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h49m26.53', '-41.30.41.1'),
          rad=17863,
          flux=   13.3,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h49m26.63', '-45.31.18.0'),
          rad=12450,
          flux=   28.5,
          bmaj=  37.4,
          bmin= 0.0,
          bpa=  46.6)

addsource(t, dir=dm.direction('J2000', '12h49m33.06', '-45.37.20.1'),
          rad=12566,
          flux=  195.5,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h49m34.40', '-48.41.07.4'),
          rad=17892,
          flux=  132.7,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h49m36.77', '-42.31.48.3'),
          rad=15523,
          flux=   21.3,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h49m49.46', '-43.10.45.7'),
          rad=14379,
          flux=   15.5,
          bmaj=  23.2,
          bmin= 0.0,
          bpa=  32.3)

addsource(t, dir=dm.direction('J2000', '12h49m50.31', '-47.39.46.6'),
          rad=15597,
          flux=   50.8,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h49m51.63', '-43.01.01.0'),
          rad=14689,
          flux=   15.8,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h49m55.80', '-42.18.47.2'),
          rad=16164,
          flux=   49.1,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h49m56.20', '-47.06.55.2'),
          rad=14579,
          flux=   17.7,
          bmaj=  25.0,
          bmin= 0.0,
          bpa= -38.3)

addsource(t, dir=dm.direction('J2000', '12h49m59.22', '-45.10.56.4'),
          rad=12706,
          flux=  325.9,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h50m03.06', '-41.41.18.5'),
          rad=17702,
          flux=   10.6,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h50m10.49', '-44.50.01.4'),
          rad=12861,
          flux=   44.3,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h50m21.94', '-45.31.31.9'),
          rad=13028,
          flux=   16.9,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h50m22.90', '-41.41.18.1'),
          rad=17862,
          flux=   24.4,
          bmaj=  48.2,
          bmin= 0.0,
          bpa= -66.7)

addsource(t, dir=dm.direction('J2000', '12h50m32.51', '-41.55.45.9'),
          rad=17358,
          flux=   24.5,
          bmaj=  19.7,
          bmin= 0.0,
          bpa= -60.9)

addsource(t, dir=dm.direction('J2000', '12h50m32.87', '-46.07.38.1'),
          rad=13556,
          flux=   24.6,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h50m34.82', '-42.06.07.3'),
          rad=16976,
          flux=   34.7,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h50m39.80', '-47.56.00.6'),
          rad=16577,
          flux=   29.7,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h50m50.39', '-43.48.14.7'),
          rad=14061,
          flux=   25.0,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h50m50.51', '-47.22.53.3'),
          rad=15540,
          flux=   29.5,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h50m52.87', '-47.52.39.7'),
          rad=16560,
          flux=   43.7,
          bmaj= 125.0,
          bmin= 0.0,
          bpa= -23.6)

addsource(t, dir=dm.direction('J2000', '12h50m53.58', '-47.52.22.8'),
          rad=16556,
          flux=   27.3,
          bmaj=  73.2,
          bmin= 0.0,
          bpa= -22.6)

addsource(t, dir=dm.direction('J2000', '12h50m54.50', '-46.35.20.2'),
          rad=14298,
          flux=  242.2,
          bmaj=  15.6,
          bmin= 0.0,
          bpa= -27.3)

addsource(t, dir=dm.direction('J2000', '12h50m57.60', '-47.50.59.4'),
          rad=16539,
          flux=   15.5,
          bmaj=  31.1,
          bmin= 0.0,
          bpa=  68.6)

addsource(t, dir=dm.direction('J2000', '12h50m59.33', '-46.32.52.8'),
          rad=14291,
          flux=   50.9,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h51m12.68', '-42.54.05.8'),
          rad=15663,
          flux=   15.5,
          bmaj=  33.1,
          bmin= 0.0,
          bpa=  12.8)

addsource(t, dir=dm.direction('J2000', '12h51m13.33', '-47.49.21.8'),
          rad=16609,
          flux=   42.5,
          bmaj=  25.5,
          bmin= 0.0,
          bpa= -60.4)

addsource(t, dir=dm.direction('J2000', '12h51m15.82', '-44.52.37.2'),
          rad=13542,
          flux=   54.1,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h51m16.99', '-48.17.56.9'),
          rad=17689,
          flux=   13.8,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h51m21.40', '-46.35.09.7'),
          rad=14551,
          flux=   12.9,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h51m23.46', '-43.50.04.0'),
          rad=14362,
          flux=   53.4,
          bmaj=  15.7,
          bmin= 0.0,
          bpa=  -2.3)

addsource(t, dir=dm.direction('J2000', '12h51m26.96', '-44.57.40.1'),
          rad=13643,
          flux=   18.3,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h51m31.11', '-43.06.54.5'),
          rad=15465,
          flux=  182.8,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h51m39.13', '-42.53.49.9'),
          rad=15921,
          flux=   73.2,
          bmaj=  14.5,
          bmin= 0.0,
          bpa=  -7.1)

addsource(t, dir=dm.direction('J2000', '12h51m43.96', '-42.29.30.1'),
          rad=16745,
          flux=   47.2,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h51m45.27', '-44.00.30.1'),
          rad=14398,
          flux=   17.0,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h51m45.68', '-46.40.43.6'),
          rad=14906,
          flux=   16.6,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h51m46.19', '-43.17.27.5'),
          rad=15328,
          flux=   13.0,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h51m46.85', '-47.49.39.4'),
          rad=16894,
          flux=   19.7,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h51m58.48', '-42.07.03.1'),
          rad=17666,
          flux=   39.0,
          bmaj=  28.2,
          bmin= 0.0,
          bpa=  83.7)

addsource(t, dir=dm.direction('J2000', '12h52m00.11', '-45.00.02.6'),
          rad=13988,
          flux=   13.0,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h52m01.07', '-42.57.43.1'),
          rad=16015,
          flux=  559.7,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h52m01.26', '-42.29.57.4'),
          rad=16887,
          flux=   30.3,
          bmaj=  52.5,
          bmin= 0.0,
          bpa=  29.1)

addsource(t, dir=dm.direction('J2000', '12h52m02.39', '-46.58.40.1'),
          rad=15493,
          flux=   83.7,
          bmaj=  15.2,
          bmin= 0.0,
          bpa=  25.5)

addsource(t, dir=dm.direction('J2000', '12h52m04.09', '-43.02.31.7'),
          rad=15906,
          flux=   74.9,
          bmaj=  31.3,
          bmin=  26.7,
          bpa=  41.0)

addsource(t, dir=dm.direction('J2000', '12h52m09.73', '-42.42.24.2'),
          rad=16561,
          flux=   12.4,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h52m14.48', '-44.42.25.3'),
          rad=14215,
          flux=   17.7,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h52m16.81', '-47.15.39.0'),
          rad=16082,
          flux=  982.9,
          bmaj=  60.2,
          bmin= 0.0,
          bpa=  67.5)

addsource(t, dir=dm.direction('J2000', '12h52m21.11', '-46.18.48.8'),
          rad=14818,
          flux=   13.0,
          bmaj=  24.2,
          bmin= 0.0,
          bpa=  71.4)

addsource(t, dir=dm.direction('J2000', '12h52m27.00', '-44.27.37.5'),
          rad=14470,
          flux=  303.2,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h52m33.15', '-47.29.15.9'),
          rad=16627,
          flux=   13.4,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h52m39.37', '-46.40.46.2'),
          rad=15420,
          flux=   13.2,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h52m47.42', '-47.54.38.6'),
          rad=17563,
          flux=   53.9,
          bmaj=  45.3,
          bmin= 0.0,
          bpa= -19.9)

addsource(t, dir=dm.direction('J2000', '12h52m47.83', '-45.00.06.8'),
          rad=14493,
          flux=   40.0,
          bmaj=  36.1,
          bmin=  29.6,
          bpa=  81.8)

addsource(t, dir=dm.direction('J2000', '12h52m48.25', '-47.40.42.8'),
          rad=17112,
          flux=  129.4,
          bmaj=  18.4,
          bmin= 0.0,
          bpa=  89.2)

addsource(t, dir=dm.direction('J2000', '12h52m54.30', '-45.02.09.1'),
          rad=14557,
          flux=   13.8,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h52m55.54', '-46.33.02.7'),
          rad=15417,
          flux=   20.9,
          bmaj=  18.8,
          bmin= 0.0,
          bpa=  60.3)

addsource(t, dir=dm.direction('J2000', '12h52m56.06', '-45.00.30.8'),
          rad=14579,
          flux=  238.7,
          bmaj=  72.3,
          bmin=  17.8,
          bpa=  82.2)

addsource(t, dir=dm.direction('J2000', '12h53m00.15', '-42.14.29.7'),
          rad=17948,
          flux=   64.4,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h53m04.50', '-44.57.35.6'),
          rad=14675,
          flux=   18.6,
          bmaj=  23.7,
          bmin= 0.0,
          bpa=  35.9)

addsource(t, dir=dm.direction('J2000', '12h53m10.75', '-45.00.54.3'),
          rad=14734,
          flux=  101.9,
          bmaj=  65.5,
          bmin= 0.0,
          bpa= -87.4)

addsource(t, dir=dm.direction('J2000', '12h53m13.74', '-46.41.43.9'),
          rad=15770,
          flux=  112.0,
          bmaj=  32.3,
          bmin= 0.0,
          bpa= -77.9)

addsource(t, dir=dm.direction('J2000', '12h53m15.76', '-45.00.20.6'),
          rad=14788,
          flux=  100.7,
          bmaj= 109.0,
          bmin=  27.7,
          bpa= -54.0)

addsource(t, dir=dm.direction('J2000', '12h53m16.03', '-43.47.01.4'),
          rad=15572,
          flux=   18.7,
          bmaj=  20.3,
          bmin= 0.0,
          bpa=  80.2)

addsource(t, dir=dm.direction('J2000', '12h53m17.47', '-45.41.47.9'),
          rad=14927,
          flux=   17.9,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h53m19.27', '-47.53.33.2'),
          rad=17791,
          flux=   16.9,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h53m20.52', '-46.33.52.5'),
          rad=15676,
          flux=   99.1,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h53m23.40', '-45.59.37.8'),
          rad=15166,
          flux=   60.8,
          bmaj=  18.1,
          bmin= 0.0,
          bpa= -63.0)

addsource(t, dir=dm.direction('J2000', '12h53m23.95', '-42.48.41.6'),
          rad=17071,
          flux=   16.8,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h53m24.62', '-43.51.28.5'),
          rad=15579,
          flux=   23.7,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h53m25.21', '-43.17.47.7'),
          rad=16298,
          flux=   43.7,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h53m37.55', '-45.05.44.4'),
          rad=15010,
          flux=   10.4,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h53m43.76', '-42.43.14.7'),
          rad=17421,
          flux=   14.1,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h53m47.70', '-45.54.51.6'),
          rad=15360,
          flux=   13.2,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h53m54.16', '-46.43.17.2'),
          rad=16190,
          flux=   14.8,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h54m01.46', '-43.32.41.4'),
          rad=16322,
          flux=   15.3,
          bmaj=  22.4,
          bmin= 0.0,
          bpa=  49.5)

addsource(t, dir=dm.direction('J2000', '12h54m07.71', '-42.57.40.0'),
          rad=17238,
          flux=   60.0,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h54m10.51', '-44.03.41.5'),
          rad=15854,
          flux=   27.3,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h54m17.36', '-43.46.04.5'),
          rad=16219,
          flux=   28.3,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h54m18.29', '-43.40.57.2'),
          rad=16326,
          flux=   18.0,
          bmaj=  23.7,
          bmin= 0.0,
          bpa= -33.0)

addsource(t, dir=dm.direction('J2000', '12h54m25.25', '-45.44.22.7'),
          rad=15650,
          flux=   15.5,
          bmaj=  25.2,
          bmin= 0.0,
          bpa=  10.4)

addsource(t, dir=dm.direction('J2000', '12h54m28.83', '-45.36.04.2'),
          rad=15629,
          flux=  194.4,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h54m34.57', '-47.29.15.4'),
          rad=17696,
          flux=   50.1,
          bmaj=  20.7,
          bmin= 0.0,
          bpa=  35.5)

addsource(t, dir=dm.direction('J2000', '12h54m37.74', '-46.31.55.5'),
          rad=16394,
          flux=   22.4,
          bmaj=  31.6,
          bmin=  22.9,
          bpa=  32.4)

addsource(t, dir=dm.direction('J2000', '12h54m38.35', '-43.25.17.0'),
          rad=16856,
          flux=   35.7,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h54m44.54', '-47.01.46.6'),
          rad=17078,
          flux=   12.8,
          bmaj=  27.4,
          bmin= 0.0,
          bpa= -36.6)

addsource(t, dir=dm.direction('J2000', '12h54m44.95', '-45.12.39.6'),
          rad=15721,
          flux=   21.0,
          bmaj=  41.9,
          bmin= 0.0,
          bpa= -31.5)

addsource(t, dir=dm.direction('J2000', '12h54m55.19', '-47.01.09.1'),
          rad=17163,
          flux=   31.4,
          bmaj=  21.9,
          bmin= 0.0,
          bpa=  46.0)

addsource(t, dir=dm.direction('J2000', '12h54m59.39', '-44.09.37.7'),
          rad=16281,
          flux=  152.5,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h55m00.34', '-43.34.50.5'),
          rad=16876,
          flux=   12.9,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h55m01.22', '-45.00.29.6'),
          rad=15902,
          flux=   13.0,
          bmaj=  24.8,
          bmin= 0.0,
          bpa=   9.8)

addsource(t, dir=dm.direction('J2000', '12h55m06.29', '-43.18.12.5'),
          rad=17298,
          flux=   29.4,
          bmaj=  18.7,
          bmin= 0.0,
          bpa=  62.3)

addsource(t, dir=dm.direction('J2000', '12h55m06.41', '-46.41.28.8'),
          rad=16852,
          flux=   76.5,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h55m07.88', '-47.22.32.4'),
          rad=17812,
          flux=   15.3,
          bmaj=  32.6,
          bmin= 0.0,
          bpa=  49.2)

addsource(t, dir=dm.direction('J2000', '12h55m08.14', '-43.34.49.5'),
          rad=16956,
          flux=   37.9,
          bmaj=  53.8,
          bmin= 0.0,
          bpa=  62.6)

addsource(t, dir=dm.direction('J2000', '12h55m09.33', '-43.38.33.6'),
          rad=16894,
          flux=   37.4,
          bmaj=  18.9,
          bmin= 0.0,
          bpa=  40.3)

addsource(t, dir=dm.direction('J2000', '12h55m40.42', '-46.29.59.8'),
          rad=16978,
          flux=   64.4,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h55m48.32', '-44.49.23.9'),
          rad=16438,
          flux=   36.2,
          bmaj=  45.7,
          bmin= 0.0,
          bpa= -21.7)

addsource(t, dir=dm.direction('J2000', '12h55m50.41', '-44.04.00.4'),
          rad=16891,
          flux=   12.2,
          bmaj=  25.1,
          bmin= 0.0,
          bpa=  87.7)

addsource(t, dir=dm.direction('J2000', '12h55m54.95', '-44.05.26.0'),
          rad=16919,
          flux=   49.1,
          bmaj=  36.0,
          bmin= 0.0,
          bpa=  63.0)

addsource(t, dir=dm.direction('J2000', '12h56m08.69', '-45.11.52.4'),
          rad=16602,
          flux=   77.6,
          bmaj=  31.5,
          bmin= 0.0,
          bpa=  36.8)

addsource(t, dir=dm.direction('J2000', '12h56m20.23', '-46.47.35.1'),
          rad=17686,
          flux=   24.6,
          bmaj=  19.3,
          bmin= 0.0,
          bpa=  17.6)

addsource(t, dir=dm.direction('J2000', '12h56m31.62', '-43.51.53.2'),
          rad=17504,
          flux=   28.3,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h56m36.74', '-46.32.29.2'),
          rad=17576,
          flux=   15.0,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h56m37.07', '-46.01.43.3'),
          rad=17165,
          flux=   41.9,
          bmaj=  17.8,
          bmin= 0.0,
          bpa=  52.1)

addsource(t, dir=dm.direction('J2000', '12h56m40.78', '-44.47.36.1'),
          rad=17002,
          flux=   21.1,
          bmaj=  46.0,
          bmin= 0.0,
          bpa=  -6.0)

addsource(t, dir=dm.direction('J2000', '12h56m42.52', '-45.57.15.7'),
          rad=17176,
          flux=   12.1,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h56m44.52', '-44.57.04.2'),
          rad=17003,
          flux=   53.8,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h57m06.99', '-44.22.60.0'),
          rad=17465,
          flux=   11.6,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h57m13.68', '-45.03.19.2'),
          rad=17295,
          flux=  362.0,
          bmaj= 0.0,
          bmin= 0.0,
          bpa=   0.0)

addsource(t, dir=dm.direction('J2000', '12h57m51.13', '-44.46.13.0'),
          rad=17752,
          flux=   80.5,
          bmaj=  20.9,
          bmin= 0.0,
          bpa= -77.4)

addsource(t, dir=dm.direction('J2000', '12h57m57.29', '-44.35.48.3'),
          rad=17884,
          flux=   33.7,
          bmaj=  18.6,
          bmin= 0.0,
          bpa=  31.6)

addsource(t, dir=dm.direction('J2000', '12h58m01.96', '-44.35.20.6'),
          rad=17936,
          flux=  929.8,
          bmaj=  33.8,
          bmin= 0.0,
          bpa= -67.6)

addsource(t, dir=dm.direction('J2000', '12h58m04.15', '-44.46.17.5'),
          rad=17889,
          flux=   47.0,
          bmaj=  15.4,
          bmin= 0.0,
          bpa= -49.8)

print t.nrows()
