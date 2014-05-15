from pylab import *

nfeeds=32

print "feeds.spacing=1.0deg"
feeds="[feed0"
for feed in range(1,4*nfeeds-1):
    feeds="%s, feed%d"%(feeds,feed)
feeds="%s, feed%d]"%(feeds,4*nfeeds-1)

print "feeds.names=%s"%feeds

feedid=0
fx=zeros(4*nfeeds, 'f')
fy=zeros(4*nfeeds, 'f')
for offx in [-0.25,0.25]:
    for offy in [-0.25,0.25]:
        for feed in range(nfeeds):
            if feed<4:
                x=feed+1
                y=0
            elif feed<28:
                x=(feed+2)%6
                y=(feed+2-x)/6
            else:
                x=feed-28+1
                y=5
            fx[feedid]=x-2.5+offx
            fy[feedid]=y-2.5+offy
            print "feeds.feed%d=[%f,%f]"%(feedid,fx[feedid],fy[feedid])
            feedid=feedid+1

