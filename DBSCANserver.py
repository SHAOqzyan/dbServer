import numpy as np
from myPYTHON import *
from mwispDBSCAN import MWISPDBSCAN



#this file used to run DBSCAN on MWISP server
#files needed

#/usr/lib/python2.7/mwispDBSCAN.py
#/usr/lib/python2.7/myPYTHON.py



doMWdbscan= MWISPDBSCAN()
doFITS=myFITS()

class dbServer(object):

    rawCOFITS=None

    def __init__(self):
        pass




    def pipeLine(self, rawCOFITS,rmsFITS=None,averageRMS=0.5,processPath="./",produceRMSFITS=False ):

        #rawCOFITS = "../fcrao_rep_m.fits"
        rawCOFITS =  rawCOFITS #"Q1Sub.fits"


        if self.rawCOFITS in None:
            return


        if produceRMSFITS and rmsFITS is None :
            rmsFITS=doFITS.getRMSFITS( self.rawCOFITS, "rms_"+self.rawCOFITS   )

        doMWdbscan.rawCOFITS =  rawCOFITS
        doMWdbscan.rmsFITS = rmsFITS #you can provde rms fits if you hae one
        doMWdbscan.averageRMS =averageRMS # if you do not have an rm fits file, use an average rms

        doMWdbscan.setDBSCANParameters( cutoff_sigma=2,minPts=4,connectivity= 1 )
        doMWdbscan.processPath =  processPath
        doMWdbscan.setCatalogSelectionCriteria( minVox=16,minChannel=3,hasBeam=1,minPeakSigma=5)


        doMWdbscan.computeDBSCAN()
        doMWdbscan.labelFITSName= doMWdbscan.getLabelFITSName()
        doMWdbscan.getCatFromLabelArray(doClean=True) #by cleaning, we remove noise clusters
        doMWdbscan.produceCleanFITS()



doServer = dbServer()


if 1:

    doServer.rawCOFITS="../MWISP_crop.fits"
    doServer.pipeLine(produceRMSFITS=True)