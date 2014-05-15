include 'note.g';
include 'gaussian.g';

# Add a gaussian
const addgaussian2d := function(a, height=1, center=[0,0], 
				fwhm=[1,1]/fwhm_to_natural, pa=0)
{

  # Rotate if necessary
  cpa := cos(pa);
  spa := sin(pa);

  width := abs(fwhm) * fwhm_to_natural;

  maxextent:=as_integer(5*fwhm[1]);
  nx:=a::shape[1];  
  ny:=a::shape[2];

  if((as_integer(center[1])<1)||(as_integer(center[1])>nx)) return a;
  if((as_integer(center[2])<1)||(as_integer(center[2])>ny)) return a;

  x:=max(1, (center[1]-maxextent)):min(nx, (center[1]+maxextent));
  y:=max(1, (center[2]-maxextent)):min(nx, (center[2]+maxextent));

  for(iy in y) {
    rx :=  cpa*(x-center[1]) + spa*(iy-center[2]);
    ry := -spa*(x-center[1]) + cpa*(iy-center[2]);
    r:=(rx/width[1])^2+(ry/width[2])^2;
    g:=height * exp(-r);
    a[x, iy, 1, ]+:=g;
  }
  return a;
}

# Add a componentlist to an image
addcomplist:=function(im, cl) {
  cs:=im.coordsys();
  ncomp:=cl.length();
  print "Processing ", ncomp, "components";
  pix:=im.getchunk();
  nx:=pix::shape[1];  
  ny:=pix::shape[2];
  dx:=180*60*60*cs.increment()[1]/pi;
  dy:=180*60*60*cs.increment()[2]/pi;
  incx:=pi*cs.increment()[1]/(180.0*60);
  incy:=pi*cs.increment()[2]/(180.0*60);
  rx:=cs.referencevalue('m', 'direction').direction.m0.value;
  ry:=cs.referencevalue('m', 'direction').direction.m1.value;
  px:=cs.referencepixel()[1];
  py:=cs.referencepixel()[2];

  for (icomp in 1:cl.length()) {
    if(icomp%100==1) print icomp;
    comp:=cl.component(icomp, T);
# The next line fails!
#    dir:=cs.topixel([comp.shape.direction.m0,
#		     comp.shape.direction.m1]);
    dir[1]:=px+(comp.shape.direction.m0.value-rx)/incx;
    dir[2]:=py+(comp.shape.direction.m1.value-ry)/incy;
    dir:=[as_integer(dir[1]), as_integer(dir[2])]
    if(is_fail(dir)) {
      print "Error converting direction ", dir::message
    }
    else if((dir[1]<1)||(dir[1]>nx)||(dir[2]<1)||(dir[2]>ny)) {
      print 'Component ', icomp, ' is off the grid: ',
	  comp.shape.direction.m0.value,
	  comp.shape.direction.m1.value, dir
    }
    else {
      if(comp.shape.type=='Point') {
	pix[dir[1], dir[2], , 1]+:=comp.flux.value;
    }
      else {
	bmaj:=dq.convert(comp.shape.majoraxis, 'arcsec').value/abs(dx);
	bmin:=dq.convert(comp.shape.minoraxis, 'arcsec').value/abs(dx);
	bpa :=dq.convert(comp.shape.positionangle, 'rad').value;
	pix:=addgaussian2d(pix, comp.flux.value[1], dir,
			   [bmaj, bmin], bpa);
      }
    }
  }
  return im.putchunk(pix);
}

include 'measures.g';
pc:=dm.direction('J2000', '12h30m00.00', '-45d00m00.0');

include 'table.g';

npix:=8192;
totalmodel := '10uJy.model';

small:=F; # Make a small image?
if(small) {
  npix:=4096;
  totalmodel := '10uJy.model.small';
}

totalcl    := '10uJy.cl';
asciifile  := 'weak.list';

tabledelete(totalmodel);
tabledelete(totalcl);
#
# Make the empty image
#  
cell:="6arcsec";
include 'coordsys.g';
cs:=coordsys(direction=T,stokes="I",spectral=T);
cs.setreferencevalue(1.4e9, 'spectral');
cs.setreferencecode('J2000', 'direction');
cs.setreferencevalue([spaste(pc.m0.value, 'rad'), spaste(pc.m1.value, 'rad')], 'direction');
cs.setincrement([cell,cell], 'direction');
cs.setreferencepixel([npix/2+1,npix/2+1], 'direction');
cs.setrestfrequency('1.420GHz');
cs.settelescope('ASKAP');

include "image.g"
im:=imagefromshape(totalmodel, shape=[npix,npix,1,1], csys=cs, overwrite=T);
im.summary();
im.done();
#
# Read the componentlist for the weak components
#
include 'componentlist.g';
cl:=emptycomponentlist();
include 'table.g';
tabledelete(cjtab);
cltab:=spaste(asciifile, '.tab');
t:=tablefromascii(cltab, asciifile);
nrows:=t.nrows();
note('Found ', nrows, ' rows');
ra:=t.getcol('RA');
dec:=t.getcol('DEC');
flux:=t.getcol('FLUX');
bmaj:=t.getcol('BMAJ');
bmin:=t.getcol('BMIN');
bpa:=t.getcol('BPA');
tflux:=0.0;
for (row in 1:nrows) {
  bpa:=bpa%360;
  cl.simulate(1, log=F);
  ncomp:=cl.length();
  cl.setflux(ncomp, [flux[row], 0.0, 0.0, 0.0], log=F);
  cl.setrefdirframe(ncomp, 'J2000', log=F);
  cl.setrefdir(ncomp, ra[row], 'deg', dec[row], 'deg', log=F);
  if((bmaj[row]>0.0)&&(bmin[row]>0.0)) {
    cl.setshape(ncomp, 'GAUSSIAN',
		majoraxis=spaste(bmaj[row], 'arcsec'),
		minoraxis=spaste(bmin[row], 'arcsec'),
		positionangle=spaste(bpa[row]%360, 'deg'), log=F);
  }
  else {
    cl.setshape(ncomp, 'POINT', log=F);
  }
  if(row%100==1) print row;
}
 
note('Found ', cl.length(), ' weak components');
include 'makenvss.g';
addnvss(cl);
note('Added strong components');
cl.rename(totalcl);

include 'image.g';
im:=image(totalmodel);
addcomplist(im, cl);
ims:=im.convolve2d(spaste(totalmodel, '.smoothed'),
		   major='30arcsec', minor='30arcsec', 
		   overwrite=T);
im.done();
ims.done();
cl.done();

exit;
