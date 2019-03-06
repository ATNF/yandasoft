fname = 'ptflux/stats.log'

f = open(fname)
try:
   pos = None
   val = None
   for line in f:
      if line.find('Max RA Dec') != -1:
         parts = line.split()
         if len(parts)<2:
            raise RuntimeError, "Expect more than one element, you have %s" % line
         val = float(parts[0])
      elif line.find('RMS MEDIAN') == -1:
         parts = line.split()
         if len(parts)!=1:
            raise RuntimeError, "Expect just one element, you have %s" % line
         pos = float(parts[0])
      else:
         if val == None or pos == None:
            raise RuntimeError, "Either argument or value cannot be parsed: val=%s, pos=%s" % (val,pos)
         print pos, val
finally:
   f.close()
