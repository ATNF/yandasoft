
# Function to compute Calculate alpha
def SimCalcAlphaBeta(imtemplate="",taylorlist=[],namealpha="",namebeta="",threshold=0.001):
   nterms = len(taylorlist);
   
   if(nterms>1):
     if(not os.path.exists(namealpha)):
       cpcmd = 'cp -r ' + imtemplate + ' ' + namealpha;
       os.system(cpcmd);
   if(nterms>2):
     if(not os.path.exists(namebeta)):
       cpcmd = 'cp -r ' + imtemplate + ' ' + namebeta;
       os.system(cpcmd);
   
   if(nterms>0):
     ia.open(taylorlist[0]);
     ptay0 = ia.getchunk();
     ia.close();
   if(nterms>1):
     ia.open(taylorlist[1]);
     ptay1 = ia.getchunk();
     ia.close();
     ia.open(namealpha);
     alpha = ia.getchunk();
     alpha.fill(0.0);
     ia.close();
   if(nterms>2):
     ia.open(taylorlist[2]);
     ptay2 = ia.getchunk();
     ia.close();
     ia.open(namebeta);
     beta = ia.getchunk();
     beta.fill(0.0);
     ia.close();

   # Calc alpha,beta from ptay0,ptay1,ptay2
   N = ptay0.shape[0];
   
   if(nterms>1):
     for ii in range(0,N):
       for jj in range(0,N):
         if(ptay0[ii,jj,0,0]>threshold):
	    mtay0 = ptay0[ii,jj,0,0];
	    mtay1 = ptay1[ii,jj,0,0];
	    alphaval = mtay1/mtay0;
	    alpha[ii,jj,0,0] = alphaval;
	    if(nterms>2):
	       mtay2 = ptay2[ii,jj,0,0];
	       beta[ii,jj,0,0] = (mtay2/mtay0) - 0.5*alphaval*(alphaval-1);
       if(ii%100 == 0):
	 print ii;

   if(nterms>1):
     ia.open(namealpha);
     ia.putchunk(alpha);
     ia.close();
   if(nterms>2):
     ia.open(namebeta);
     ia.putchunk(beta);
     ia.close();



##############################################################################


msname = 'ptest.ms';
rname = 'try';
niter=100;
scales=[0];
ntaylor=3;
threshold='0.05mJy';
imsize=1024;
cellsize="8.0arcsec";
stokes="I";
reffreq="1.4GHz";
gain=0.5;
mask="";
algo="msmfs";
weighting="briggs";
ftm="ft";


models=[];
restoreds=[];
residuals=[];
masks=[];

im.close();
im.open(msname);
im.selectvis(spw=0);
if(algo=="msmfs"):
   for tt in range(0,ntaylor):
        models.append(rname+'.'+str(tt)+'.model');
        restoreds.append(rname+'.'+str(tt)+'.restored');
        residuals.append(rname+'.'+str(tt)+'.residual');
        masks.append(mask);
        im.defineimage(nx=imsize,ny=imsize,cellx=cellsize,celly=cellsize,nchan=1,stokes=stokes,mode='mfs');
        im.make(models[tt]);
else:
   models = [rname+'.model'];
   restoreds = [rname+'.restored'];
   residuals = [rname+'.residual'];
   masks = mask;
   im.defineimage(nx=imsize,ny=imsize,cellx=cellsize,celly=cellsize,nchan=1,stokes=stokes,mode='mfs');
   im.make(models[0]);

im.weight(type=weighting);
im.setscales(scalemethod='uservector',uservector=scales);
im.settaylorterms(ntaylorterms=ntaylor,reffreq=(qa.convert(qa.unit(reffreq),"Hz"))['value']);
im.setoptions(ftmachine=ftm);
if(mask == ""):
  print 'clean without mask';
  im.clean(model=models,image=restoreds,residual=residuals,algorithm=algo,threshold=threshold,niter=niter,interactive=False,gain=gain); 
else:
  print 'clean with mask';
  im.clean(model=models,image=restoreds,residual=residuals,algorithm=algo,threshold=threshold,niter=niter,interactive=False,gain=gain,mask=masks); 

im.done();
###########################################################################

## Calculate alpha and beta

imtemplate=rname+'.0.restored';
taylist=[];
for i in range(0,ntaylor):
	taylist.append(rname+'.'+str(i)+'.restored');

SimCalcAlphaBeta(imtemplate=imtemplate,taylorlist=taylist,namealpha=rname+'.'+ftm+'.restored.alpha',namebeta=rname+'.'+ftm+'.restored.beta',threshold=0.02);


###########################################################################

