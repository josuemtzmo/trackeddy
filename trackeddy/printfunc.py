## All this functions were grenerated to dynamically print output.
import sys
import numpy as np
import time

class Printer():
    """
    Print things to stdout on one line dynamically
    """ 
    def __init__(self):
        self.tic=time.clock()
        sys.stdout.flush()
        
    def printtextoneline(self,string):
        sys.stdout.write("\r\x1b[K"+string.__str__())
        sys.stdout.flush()
        
    def timepercentprint(self,minv,maxv,step,i,neddies):
        percent=(float(i+1)/maxv)*100.0
        etime=round(time.clock()-self.tic)
        stmtime=round((etime/percent)*100)
        
        progress=int(20/(maxv/(step*(i+1))))
        emptyprog=20-progress
        sys.stdout.write("\r 0% [{0}>{1}]{2}% | Elapsed Time:{3} s | Estimated Time:{4} s | Info: {5} |".format("="*progress,' '*emptyprog,round(percent),etime,stmtime,neddies))
        if percent != 100:
            sys.stdout.flush()
        else:
            print('')
          
