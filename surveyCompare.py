from myPYTHON import *

import numpy as np

from mwispDBSCAN import MWISPDBSCAN



doFITS = myFITS()


class surveyCom(object):

    dataPath="/media/qzyan/maclinux/Data/Q2CompareData/"

    #lRange= [ 104.75, 141.54  ]
    lRange= [ 110, 141.54  ]

    #bRange =[ -3.028, 5.007 ]
    bRange =[ -2, 4 ]

    vRange = [-137, 32  ]


    surveyNames = ["MWISP","CfA","OGS"]

    surveyRawFITS=[dataPath+"100_150_U.fits", dataPath+"cfa_rep.fits", dataPath+"fcrao_reSave.fits"  ]
    surveyCropFITS=[dataPath+"MWISP_crop.fits", dataPath+"CfA_crop.fits", dataPath+"OGS_crop.fits"  ]
    surveyCropRMSFITS=[ "MWISP_crop_RMS.fits",  "CfA_crop_RMS.fits", "OGS_crop_RMS.fits"  ]



    surveyN=3

    def __init__(self):
        pass



    def prepareData(self):
        """
        crop the data, then run DBSCAN
        :return:
        """


        for i in range(self.surveyN):

            rawFITS=self.surveyRawFITS[i]
            outFITS=self.surveyCropFITS[i]
            doFITS.cropFITS(rawFITS,outFITS=outFITS,Vrange=self.vRange,Lrange=self.lRange,Brange=self.bRange,overWrite=True )




    def getRMSFITS(self):
        """
        calculate the RMS files of the four
        :return:
        """
        for i in range(self.surveyN):

            cropFITS = self.surveyCropFITS[i]
            saveRMSFITS = self.surveyCropRMSFITS[i]

            doFITS.getRMSFITS(cropFITS, saveRMSFITS)



    def findClouds(self):
        """
        fidn clouds
        :return:
        """

        for i in Range(self.surveyN):
            print i

        doMWdbscan = MWISPDBSCAN()
        doMWdbscan.rawCOFITS = None
        doMWdbscan.rmsFITS = self.getRMSFITS()

        doMWdbscan.setDBSCANParameters( cutoff_sigma=cutoff,minPts=minPts,connectivity=contype)
        doMWdbscan.processPath = self.tmpPath

        doMWdbscan.computeDBSCAN()
        doMWdbscan.getCatFromLabelArray(doClean=True)

        #only produce clean fits for smFactor 1.0
        #if smFactor==1.0:
        doMWdbscan.produceCleanFITS()


doSurvey = surveyCom()

#doSurvey.prepareData()

doSurvey.getRMSFITS()