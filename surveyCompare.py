from myPYTHON import *
import matplotlib.pyplot as plt

import numpy as np
import matplotlib as mpl
from mwispDBSCAN import MWISPDBSCAN



doFITS = myFITS()


class surveyCom(object):

    dataPath="/media/qzyan/maclinux/Data/Q2CompareData/"
    DBSCANresults="./DBSCANresults/"

    #lRange= [ 104.75, 141.54  ]
    lRange= [ 110, 141.54  ]

    #bRange =[ -3.028, 5.007 ]
    bRange =[ -2, 4 ]

    vRange = [-137, 32  ]


    surveyNames = ["MWISP","CfA","OGS"]

    Imwisp=0
    ICfA=1
    IOGS=2

    surveyRawFITS=[dataPath+"100_150_U.fits", dataPath+"cfa_rep.fits", dataPath+"fcrao_reSave.fits"  ]
    surveyCropFITS=[dataPath+"MWISP_crop.fits", dataPath+"CfA_crop.fits", dataPath+"OGS_crop.fits"  ]
    surveyCropRMSFITS=[ "MWISP_crop_RMS.fits",  "CfA_crop_RMS.fits", "OGS_crop_RMS.fits"  ]

    tbSuffix="_cropdbscanS2P4Con1_Clean.fit"
    surveyTableNames= [ DBSCANresults+surveyNames[0] +  tbSuffix , DBSCANresults+surveyNames[1] +  tbSuffix , DBSCANresults+surveyNames[2] +  tbSuffix ]

    surveyTBs=[]

    for eachTBName in surveyTableNames:
        surveyTBs.append(Table.read(eachTBName))

    surveyN=3

    surveyPixelResolution = [30./60,0.125*60, 0.01395*60 ] #in arcmin
    surveyVelResolution = [ 0.167,    1.300, 0.812  ] #in km/s


    def __init__(self):
        pass



    def prepareData(self):
        """
        crop the data, then run DBSCAN
        :return:
        """


        for i in range(self.surveyN):

            if i ==0:
                continue
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



    def checkTable(self):
        """

        :return:
        """

        for eachTB in self.surveyTBs:

            print len(eachTB)


    def getBinAndCenter(self,a,binsEdges ):
        """
        here, bins need to an array
        :param a:
        :param bins:
        :return:
        """
        histArray, bin_edges =  np.histogram(a,binsEdges )

        binCenters = ( bin_edges[0:-1] + bin_edges[1:] )*0.5
        histArray= histArray /1.0
        #histArray[histArray==0]=np.nan


        return histArray,binCenters

    def drawLhistgram(self,binN=20,colName="x_cen",colLabel="Galactic longitude (deg)",saveLabel=""):

        """
        drawThe histgram of molecular clouds along the Galactic longitude
        :return:
        """

        if colName=="x_cen":

            bins=np.linspace( np.min(self.lRange) ,np.max(self.lRange), binN  )
        else:
            bins=np.linspace( np.min(self.bRange) ,np.max(self.bRange), binN  )


        #n, bins, patches = plt.hist(x, num_bins, facecolor='blue', alpha=0.5)
        lCol= colName


        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot(1, 1, 1)
        # fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(12,6),sharex=True)
        rc('text', usetex=True)
        rc('font', **{'family': 'sans-serif', 'size': 19, 'serif': ['Helvetica']})
        mpl.rcParams['text.latex.preamble'] = [
        r'\usepackage{tgheros}', # helvetica font
        r'\usepackage{sansmath}', # math-font matching helvetica
        r'\sansmath' # actually tell tex to use it!
        r'\usepackage{siunitx}', # micro symbols
        r'\sisetup{detect-all}', # force siunitx to use the fonts
        ]
        #ax.scatter( onSourceIndex[1], onSourceIndex[0],s=5,color="gray")


        histArray,centerS = self.getBinAndCenter(self.surveyTBs[0][ lCol ],bins)
        ax.plot(centerS, histArray, 'o-', color='blue')


        #ax.hist(self.surveyTBs[0][ lCol ] ,bins,   facecolor="blue" ,zorder=1 ,alpha=0.5,label=self.surveyNames[0])
        #ax.hist(self.surveyTBs[2][ lCol ] ,bins,   facecolor="green" ,zorder=2 ,alpha=0.5 ,label=self.surveyNames[2])
        #ax.hist(self.surveyTBs[1][ lCol ] ,bins,   facecolor="red" ,zorder=3 ,alpha=0.5 ,label=self.surveyNames[1])


        ax.legend(loc=1)

        ax.set_xlabel( colLabel )

        ax.set_ylabel(r"Number of molecular clouds")


        plt.savefig(saveLabel+"histNumber.png", bbox_inches='tight',dpi=600)

    def selectTBbyLBVRange(self,TB,lRange,bRange,vRange=None):



        selectCriteria = TB["x_cen"]>=min(lRange)

        selectCriteria=np.logical_and(  selectCriteria,TB["x_cen"]<=max(lRange)    )

        selectCriteria = np.logical_and(selectCriteria, TB["y_cen"] <= max(bRange))
        selectCriteria = np.logical_and(selectCriteria, TB["y_cen"] >= min(bRange))

        if vRange is not None:
            selectCriteria = np.logical_and(selectCriteria, TB["v_cen"] <= max(vRange))
            selectCriteria = np.logical_and(selectCriteria, TB["v_cen"] >= min(vRange))


        return TB[selectCriteria]

    def gridCompare(self, lN=10,bN=4,vN=4):



        lEdges = np.linspace( np.min(self.lRange) ,np.max(self.lRange), lN+1  )

        bEdges = np.linspace( np.min(self.bRange) ,np.max(self.bRange), bN+1  )
        vEdges = np.linspace( np.min(self.vRange) ,np.max(self.vRange), vN+1  )




        mwispN = []
        tbCfAN = []
        tbOGSN = []
        lRanges=[]
        bRanges = []
        vRanges = []

        for i in range(lN-1):
            for j in range(bN-1):
                for k in range(vN - 1):

                    countLrange= lEdges[i:i+2]
                    countBrange= bEdges[j:j+2]
                    cloutVrange = vEdges[k:k+2]

                    tbMWISP = self.selectTBbyLBVRange(self.surveyTBs[self.Imwisp],countLrange,countBrange,vRange=cloutVrange)
                    tbCfA = self.selectTBbyLBVRange(self.surveyTBs[self.ICfA],countLrange,countBrange,vRange=cloutVrange)
                    tbOGS = self.selectTBbyLBVRange(self.surveyTBs[self.IOGS ],countLrange,countBrange,vRange=cloutVrange)

                    mwispN.append(  len(tbMWISP)/1. )
                    tbCfAN.append(  len(tbCfA)/1. )
                    tbOGSN.append(  len(tbOGS)/1. )

                    lRanges.append(  countLrange )

                    bRanges.append(  countBrange )
                    vRanges.append( cloutVrange )
        mwispN = np.asarray(mwispN)
        tbCfAN = np.asarray(tbCfAN)
        tbOGSN = np.asarray(tbOGSN)

        #####################

        numberN = np.arange(len(mwispN))

        ratio= mwispN/tbOGSN
        ratio[tbOGSN==0] = np.nan
        maxIndex=np.nanargmax(ratio)


        #maxIndex = mwispN.argmax()

        print lRanges[maxIndex]
        print bRanges[maxIndex]
        print vRanges[maxIndex]

        #cropfits

        for i in range(self.surveyN):

            eachRawFITS  = self.surveyCropFITS[i]
            doFITS.cropFITS(eachRawFITS, Lrange=lRanges[maxIndex], Brange=bRanges[maxIndex], Vrange =  vRanges[maxIndex] ,outFITS= self.surveyNames[i]+"MaximumCompare.fits",overWrite=True   )





        print mwispN[maxIndex]
        print tbCfAN[maxIndex]
        print tbOGSN[maxIndex]

        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot(1, 1, 1)
        # fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(12,6),sharex=True)
        rc('text', usetex=True)
        rc('font', **{'family': 'sans-serif', 'size': 19, 'serif': ['Helvetica']})
        mpl.rcParams['text.latex.preamble'] = [
        r'\usepackage{tgheros}', # helvetica font
        r'\usepackage{sansmath}', # math-font matching helvetica
        r'\sansmath' # actually tell tex to use it!
        r'\usepackage{siunitx}', # micro symbols
        r'\sisetup{detect-all}', # force siunitx to use the fonts
        ]
        #ax.scatter( onSourceIndex[1], onSourceIndex[0],s=5,color="gray")
        ax.step( numberN,mwispN,where="mid" ,color="blue",lw=1 ,label=self.surveyNames[self.Imwisp])
        ax.step( numberN,tbCfAN,where="mid" ,color="red",lw=1  ,label=self.surveyNames[self.ICfA])
        ax.step( numberN,tbOGSN,where="mid" ,color="green",lw=1 ,label=self.surveyNames[self.IOGS])

        ax.legend(loc=1,handlelength=0.8)

        ax.set_ylabel(r"Number of molecular clouds")
        ax.set_xlabel(r"Region index")

        ax.axvline(x=maxIndex , ls="-", color='black', lw=1.5)
        plt.savefig(  "histNumber_grid.png", bbox_inches='tight',dpi=600)


    def getFlux(self,surveyID):
        """

        :param surveyID:
        :return:
        """

        return self.surveyTBs[surveyID]["sum"]*self.surveyVelResolution[surveyID]*self.surveyPixelResolution[surveyID]**2



    def compareFlux(self,binN=1000):
        """

        :param binStep:
        :return:
        """



        areaMWISP =  self.getFlux(self.Imwisp )
        areaCfA  =   self.getFlux(self.ICfA ) #self.surveyTBs[self.ICfA]["sum"]
        areaOGS =   self.getFlux(self.IOGS ) #self.surveyTBs[self.IOGS]["sum"]


        minArea= min( [ np.min( areaMWISP )   ,  np.min( areaCfA )    ,  np.min( areaOGS )    ] )
        maxArea= max( [ np.max( areaMWISP )   ,  np.max( areaCfA )    ,  np.max( areaOGS )    ] )




        #areaBinEdges = np.linspace(minArea ,maxArea ,  binN  )
        areaBinEdges = np.linspace(1 ,10000 ,  binN  )

        #Compare the histgram of molecular clouds along area
        fontsize= 18
        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot(1, 1, 1)
        # fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(12,6),sharex=True)
        rc('text', usetex=True)
        rc('font', **{'family': 'sans-serif', 'size': fontsize , 'serif': ['Helvetica']})
        mpl.rcParams['text.latex.preamble'] = [
        r'\usepackage{tgheros}', # helvetica font
        r'\usepackage{sansmath}', # math-font matching helvetica
        r'\sansmath' # actually tell tex to use it!
        r'\usepackage{siunitx}', # micro symbols
        r'\sisetup{detect-all}', # force siunitx to use the fonts
        ]
        #ax.scatter( onSourceIndex[1], onSourceIndex[0],s=5,color="gray")


        histArray,centerS = self.getBinAndCenter( areaMWISP ,areaBinEdges)
        ax.plot(centerS, histArray, 'o-', color='blue',lw=1,markersize=3,label=self.surveyNames[self.Imwisp])


        histArray,centerS = self.getBinAndCenter( areaOGS ,areaBinEdges)
        ax.plot(centerS, histArray, 'o-', color='green',lw=1,markersize=3,label=self.surveyNames[self.IOGS])

        histArray,centerS = self.getBinAndCenter( areaCfA ,areaBinEdges)
        ax.plot(centerS, histArray, 'o-', color='red',lw=1,markersize=3,label=self.surveyNames[self.ICfA])


        #ax.hist( areaMWISP ,areaBinEdges,   facecolor="blue" ,zorder=1 ,alpha=0.8,label=self.surveyNames[self.Imwisp])
        #ax.hist( areaOGS ,areaBinEdges,   facecolor="green" ,zorder=2 ,alpha=0.8,label=self.surveyNames[self.IOGS])
        #ax.hist( areaCfA ,areaBinEdges,   facecolor="red" ,zorder=3 ,alpha=0.8,label=self.surveyNames[self.ICfA])

        ax.tick_params(axis='both', which='major', labelsize=fontsize)
        ax.tick_params(axis='both', which='minor', labelsize=fontsize)

        ax.set_yscale('log')

        ax.set_xscale('log')
        #ax.set_xlim([0.5,10000])
        ax.legend(loc=1)
        ax.set_ylabel(r"Number of molecular clouds",fontsize=fontsize)
        ax.set_xlabel(r"Flux (K km/s arcmins)",fontsize=fontsize)



        plt.savefig(  "hist_flux.png", bbox_inches='tight',dpi=600)








    def compareArea(self,binN=1000):
        """

        :param binStep:
        :return:
        """
        areaMWISP = self.surveyTBs[self.Imwisp]["area_exact"]
        areaCfA  = self.surveyTBs[self.ICfA]["area_exact"]
        areaOGS = self.surveyTBs[self.IOGS]["area_exact"]


        minArea= min( [ np.min( areaMWISP )   ,  np.min( areaCfA )    ,  np.min( areaOGS )    ] )
        maxArea= max( [ np.max( areaMWISP )   ,  np.max( areaCfA )    ,  np.max( areaOGS )    ] )


        #areaBinEdges = np.linspace(minArea ,maxArea ,  binN  )
        areaBinEdges = np.linspace(1 ,10000 ,  binN  )

        #Compare the histgram of molecular clouds along area
        fontsize= 18
        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot(1, 1, 1)
        # fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(12,6),sharex=True)
        rc('text', usetex=True)
        rc('font', **{'family': 'sans-serif', 'size': fontsize , 'serif': ['Helvetica']})
        mpl.rcParams['text.latex.preamble'] = [
        r'\usepackage{tgheros}', # helvetica font
        r'\usepackage{sansmath}', # math-font matching helvetica
        r'\sansmath' # actually tell tex to use it!
        r'\usepackage{siunitx}', # micro symbols
        r'\sisetup{detect-all}', # force siunitx to use the fonts
        ]
        #ax.scatter( onSourceIndex[1], onSourceIndex[0],s=5,color="gray")




        histArray,centerS = self.getBinAndCenter( areaMWISP ,areaBinEdges)
        ax.plot(centerS, histArray, 'o-', color='blue',lw=1,markersize=3,label=self.surveyNames[self.Imwisp])


        histArray,centerS = self.getBinAndCenter( areaOGS ,areaBinEdges)
        ax.plot(centerS, histArray, 'o-', color='green',lw=1,markersize=3,label=self.surveyNames[self.IOGS])

        histArray,centerS = self.getBinAndCenter( areaCfA ,areaBinEdges)
        ax.plot(centerS, histArray, 'o-', color='red',lw=1,markersize=3,label=self.surveyNames[self.ICfA])

        #ax.hist( areaMWISP ,areaBinEdges,   facecolor="blue" ,zorder=1 ,alpha=0.8,label=self.surveyNames[self.Imwisp])
        #ax.hist( areaOGS ,areaBinEdges,   facecolor="green" ,zorder=2 ,alpha=0.8,label=self.surveyNames[self.IOGS])
        #ax.hist( areaCfA ,areaBinEdges,   facecolor="red" ,zorder=3 ,alpha=0.8,label=self.surveyNames[self.ICfA])






        ax.tick_params(axis='both', which='major', labelsize=fontsize)
        ax.tick_params(axis='both', which='minor', labelsize=fontsize)

        ax.set_yscale('log')

        ax.set_xscale('log')
        ax.set_xlim([0.5,10000])
        ax.legend(loc=1,handlelength=0.8)
        ax.set_ylabel(r"Number of molecular clouds",fontsize=fontsize)
        ax.set_xlabel(r"Angular area (arcmins)",fontsize=fontsize)



        plt.savefig(  "hist_area.png", bbox_inches='tight',dpi=600)


    def compareRatios(self):
        """
        compare the ratio of
        :return:
        """



        areaMWISP = self.surveyTBs[self.Imwisp]["area_exact"]
        areaCfA  = self.surveyTBs[self.ICfA]["area_exact"]
        areaOGS = self.surveyTBs[self.IOGS]["area_exact"]

        fluxMWISP =  self.getFlux(self.Imwisp )
        fluxCfA  =   self.getFlux(self.ICfA ) #self.surveyTBs[self.ICfA]["sum"]
        fluxOGS =   self.getFlux(self.IOGS ) #self.surveyTBs[self.IOGS]["sum"]



        areaMWISP,areaOGS, areaCfA  =    np.sum(areaMWISP), np.sum(areaOGS) , np.sum(areaCfA)
        fluxMWISP, fluxOGS,  fluxCfA =  np.sum(fluxMWISP) , np.sum(fluxOGS) , np.sum(fluxCfA)


        print "MWISP/OGS", "MWISP/CfA"

        print "Area: {:.1f}    {:.1f}".format(areaMWISP/areaOGS,areaMWISP/areaCfA   )
        print "Flux: {:.1f}    {:.1f}".format(fluxMWISP/fluxOGS,fluxMWISP/fluxCfA   )



    def ZZZ(self):
        """

        :return:
        """

doSurvey = surveyCom()


doSurvey.compareRatios()
#doSurvey.compareFlux()
#doSurvey.compareArea()

#doSurvey.gridCompare()


if 0:
    #doSurvey.drawLhistgram(saveLabel="L_")
    #doSurvey.drawLhistgram(saveLabel="B_",colName="y_cen",colLabel="Galactic latitude (deg)" ,binN=30  )

    doSurvey.prepareData()

    #doSurvey.getRMSFITS()