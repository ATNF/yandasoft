#
# Define xNTD configuration for use in AIPS++ simulator
#
readconf:=function(filename='xntd.conf') {

  include 'measures.g';

  rec:=[=];

  rec.location:=dm.position('WGS84', '+117.471deg', '-25.692deg');
  
  f:=open(spaste('< ', filename));
  line:=read(f);
  lineno:=0;
  nant:=0;
  while(len(line)) {
    parts:=split(line, ',');
    if(len(parts)==3) {
      nant+:=1;
      rec.x[nant]:=as_float(parts[1]);
      rec.y[nant]:=-1.0*as_float(parts[2]);
      rec.z[nant]:=0.0;
    }
    lineno+:=1;
    line:=read(f);
  }
  f:=F;
  rec.mount:=array('equatorial', nant);
  rec.diam:=array(12.0, nant);
  for (i in 1:nant) {
    rec.name[i]:=sprintf('MIRANDA%02d', i);
  }
  return rec;
}

plotconf:=function(filename,baseline=3500) {
  include 'pgplotter.g';
  p:=pgplotter(background='white', foreground='black');
  p.subp(2,1);
  r:=readconf(filename);
  r.x-:=sum(r.x)/len(r.x);
  r.y-:=sum(r.y)/len(r.y);
  nant:=len(r.x);
  p.env(-baseline/2,baseline/2,-baseline/2,baseline/2,1,1);
  p.setplottitle(spaste('xNTD configuration: ', filename));
  p.setxaxislabel('X (m)');
  p.setyaxislabel('Y (m)');
  p.iden();
  p.sci(2);
  p.pt(r.x, r.y, 4);
  p.sci(1);
  p.env(-baseline,baseline,-baseline,baseline,1,1);
  p.setplottitle(spaste('xNTD configuration: ', filename));
  p.setxaxislabel('U (m)');
  p.setyaxislabel('V (m)');
  p.iden();
  p.sci(2);
  u:=array(0.0, nant, nant);
  v:=array(0.0, nant, nant);
  for (i in 1:nant) {
    u[,i]:=r.x;
    v[,i]:=r.y;
  }
  for (i in 1:nant) {
    u[i,]-:=r.x;
    v[i,]-:=r.y;
  }
  p.pt(u, v, 2);
  p.sci(1);
  return p;
}

