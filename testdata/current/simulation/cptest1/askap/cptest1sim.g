#
# AIPS++ script to simulate ASKAP observations at 1.5GHz
#
#
# Parameters
#
nx:=1; # Number of feeds in x
ny:=1; # Number of feeds in y
config:='askap30';
nchan:=128; # total number of channels in 256MHz
nspw:=4 # Number of spectral windows over which the channels are spread
#
npix:=8192;
totalmodel := '10uJy.model';

small:=F; # Make a small image?
if(small) {
  npix:=4096;
  totalmodel := '10uJy.model.small';
}


include 'logger.g';
dl.screen();
include 'note.g';
include 'vpmanager.g';
include 'image.g';
  

dl.purge(0);


a:=array(0.0, 2*nx*ny);
if(nx*ny>1) {
  i:=1;
  for (y in 1:ny) {
    for (x in 1:nx) {
      a[2*i-1]:=as_float(x)-as_float(nx+1)/2.0;
      a[2*i]:=as_float(y)-as_float(ny+1)/2.0;
      i+:=1;
    }
  }
}
else {
  a[1]:=0.0;
  a[2]:=0.0;
}
print a;
nfeeds:=len(a)/2;
#
# Get the array configuration
#
include 'readconf.g';
rec:=readconf(spaste(config, '.txt'));
print rec;
nant:=len(rec.x);
rec.mount:=array('equatorial', nant);
#
# Spacing between pointings in radians
#
spacing:=pi*30/(60*180);
#
cell:='6arcsec';

#
# Integration time, time on source, number of integrations
#
t:=30; # Integration time
tsource:=8*3600; # Total time on source
thourangle:=8*3600.0; # Range in hour angle
tstep:=30; # Integrate for this time
#nint:=as_integer(0.5+tstep/t); # Number of integrations per step
nint:=1; # Number of integrations per step
nstep:=as_integer(0.5+tsource/tstep); # Number of steps
tjump:=thourangle/nstep;

integrationtime:=spaste(t,'s');

include 'measures.g';
pc:=dm.direction('J2000', '12h30m00.00', '-45d00m00.0');

deltafreq:=spaste(256.0/nchan,'MHz');

note('***** xNTD Observing details *****');
note('Integration time       = ', integrationtime);
note('Number of integrations = ', nint);
note('Number of steps        = ', nstep);
note('Number of channels     = ', nchan);
note('Number of spw          = ', nspw);
note('Channel width          = ', deltafreq);
note('Cellsize               = ', cell);
  
for (spwid in 1:nspw) {

  testdir := spaste('data/spw_', spwid);
  msname     := spaste(testdir, '/sim.ms');
  simvp      := spaste(testdir, '/sim.vp');
#
# Do we want to simulate new data?
#
  shell(paste('rm -rf ', testdir));
  ok := shell(paste("mkdir -p", testdir));
  if (ok::status) { throw("mkdir", testdir, "fails!") };
#
# Need a voltage pattern description. This is an Airy disk
# for a 12m diameter antenna with uniform illumination and
# blockage of 1m diameter.
#
  vp:=vpmanager();
  print vp.setpbairy('ASKAP', dopb=T, dishdiam='12m', blockagediam='1m', maxrad='10deg');
  vp.saveastable(simvp);
  
  note('Create the empty measurementset');
  include "newsimulator.g";
#
# CONRAD measures tables are set for this date
#    
  reftime := dm.epoch('utc', '8Jun2007');
#
# Define a (new)simulator tool
#
  mysim := newsimulator(msname);
#
# Set options
#
#
# Set up the configuration
#
  print mysim.setconfig('ASKAP', x=rec.x[1:nant], y=rec.y[1:nant], z=rec.z[1:nant], 
			dishdiameter=rec.diam[1:nant], mount=rec.mount[1:nant],
			antname=rec.name[1:nant], coordsystem='local',
			referencelocation=rec.location);

#
  for (feed in 1:nfeeds) {
    feedposx:=spaste(pc.m0.value+a[2*feed-1]*spacing/cos(pc.m1.value), 'rad');
    feedposy:=spaste(pc.m1.value+a[2*feed]*spacing,   'rad');
    mysim.setfield(spaste('ASKAP_FEED_', feed),
		   sourcedirection=dm.direction('j2000', feedposx, feedposy));
  }
#
# Spectral window: 
#
  dfreq:=0.256/as_float(nspw);
  freq0:=1.420-0.256+as_float(spwid-1)*dfreq;
  mysim.setspwindow(spwname=spaste('LBAND', spwid), 
		    freq=spaste(freq0, 'GHz'),
		    deltafreq=deltafreq,
		    freqresolution=deltafreq,
		    nchannels=nchan/nspw,
		    stokes='XX XY YX YY');
#
# Various limits, etc.
#    
  mysim.setlimits(shadowlimit=0.001, elevationlimit='8.0deg');
  mysim.setauto(autocorrwt=0.0);
  mysim.setfeed('perfect X Y', pol="");
  
  mysim.settimes(integrationtime=integrationtime, usehourangle=T,
		 referencetime=reftime);
#
# Loop over feed. 
#    
  tstart:=-thourangle/2.0;
  for (step in 1:nstep) {
    for (feed in 1:nfeeds) {
      print "Observing ", spaste('ASKAP_FEED_', feed), spaste('LBAND',spwid), spaste(tstart, 's'), spaste(tstart+nint*t, 's');
      if(!mysim.observe(spaste('ASKAP_FEED_', feed), 
			spaste('LBAND',spwid),
			starttime=spaste(tstart, 's'),
			stoptime=spaste(tstart+nint*t, 's'))) {
	return throw('observe failed');
      }
    }
    tstart+:=tjump;
  }
#
# Add noise. We add the model using imager since simulator seems to be
# broken
#
  print mysim.summary();

  dl.printtofile(filename='cptest1sim.log');
  mysim.setnoise('calculate', antefficiency=0.65, correfficiency=0.99,
		 spillefficiency=0.95, tau=0.0, trx=50.0,
		 tatmos=250.0, tcmb=2.7);
  print mysim.setseed(random());
  mysim.corrupt();
  mysim.done();
  include 'imager.g';
  myim:=imager(msname);
  myim.setimage();
  print myim.setoptions(ftmachine='wproject', wprojplanes=64);
  print myim.setvp(dovp=T, usedefaultvp=F, vptable=simvp, parangleinc="360deg");
  print myim.ft(totalmodel);
  myim.done();
  ot:=table(msname, readonly=F);
  md:=ot.getcol('MODEL_DATA');
  print max(abs(md));
  md+:=ot.getcol('DATA');
  print ot.putcol('DATA', md);
  print ot.putcol('CORRECTED_DATA', md);
  ot.close();
  dl.printtofile(filename='cptest1sim.log');
}

exit;

