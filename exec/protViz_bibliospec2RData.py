#!/usr/bin/python
# -*- coding: latin1 -*-


# Christian Panse <cp@fgcz.ethz.ch>

# Copyright (C) 2014 Functional Genomics Center Zurich ETHZ|UZH. All rights reserved.
# 
# # Licensed under  GPL version 3
# $HeadURL: http://fgcz-svn.unizh.ch/repos/fgcz/testing/proteomics/R/protViz/exec/protViz_bibliospec2RData.py $
# $Id: protViz_bibliospec2RData.py 6318 2014-04-02 13:02:51Z cpanse $


"""
# INPUT: bibliospec sqlite files
# OUTPUT: extracts MS2 
"""

import sqlite3
import zlib
import sys
import numpy
import os
import math
import re
import string



def myquery(sqllitedb):
    count=0
    myquery = "SELECT numPeaks, peakMZ, peakIntensity, peptideSeq, precursorCharge, precursorMZ, retentionTime, peptideModSeq, score, SpectrumSourceFiles.fileName FROM SpectrumSourceFiles, RefSpectraPeaks, RefSpectra WHERE RefSpectra.id=RefSpectraPeaks.RefSpectraID and SpectrumSourceFiles.id = RefSpectra.fileID;"
    myconnect = sqlite3.connect(sqlitedb)
    mycursor = myconnect.cursor()

    myRDataName = os.path.basename(sqllitedb)
    myRFile = open(sqllitedb+".R", 'w')
    myRFile.write(myRDataName + "<-list()")

    for row in mycursor.execute(myquery):
        count = count + 1

        # print count
        try:
            # print " COMPRESSED ********************"
            myblob = zlib.decompress(row[1]) 
            MZ = numpy.fromstring(myblob, dtype='d')
        except zlib.error as e:
            # print " NOT COMPRESSED ********************"
            # print (e)
            MZ = numpy.fromstring(row[1], dtype='d')
            # print "[DONE]"
        except:
            print "+++ unexpected error +++"
            sys.exit(1)

        try:
            # print " COMPRESSED ********************"
            myblob = zlib.decompress(row[2]) 
            Intensity=numpy.fromstring(myblob, dtype='f')

        except zlib.error as e:
            # print " NOT COMPRESSED ********************"
            print (e)
            Intensity=numpy.fromstring(row[2], dtype='f')
        except:
            print "+++ unexpected error +++"
            sys.exit(1)


        # print "# SCAN =", count, "|", row[3], "|", row[4], "|", row[5],"|", row[6]
        # for idx in range(1, len(MZ)-1):
        #    print (MZ[idx]), ",", 

        myRFile.write( "\n\n\n" + myRDataName + "[[" + '{}'.format(count) + "]] <- list(" )

        myRFile.write("\n\tmZ=c(")
        for idx in range(len(MZ)):
            myRFile.write('{}'.format(MZ[idx]))
            if idx < len(MZ)-1:
                myRFile.write(", ")
        myRFile.write("),\n")

        myRFile.write("\n\tintensity=c(")
        for idx in range(len(MZ)):
            myRFile.write('{}'.format(Intensity[idx]))
            if idx < len(MZ)-1:
                myRFile.write(", ")
        myRFile.write("),\n")

        myRFile.write("\tpeptideSequence='" + str(row[3]) + "',\n")
        myRFile.write("\tcharge=" + str(row[4]) + ",\n")
        myRFile.write("\tpepmass=" + str(row[5]) + ",\n")
        # there is no mascot score defined. it looks more like an E-value

        myRFile.write("\tpeptideModSeq='"  + str(row[7]) + "',\n")
        varModification = row[7]
        varModification = re.sub("[A-Z]", "0,", varModification)
        varModification = re.sub("0,\[", "", varModification)
        varModification = re.sub("\]", ",", varModification)
        varModification = re.sub(",$", "", varModification)
        varModification = "c(" + varModification  + ")"

        myRFile.write("\tvarModification="  + str(varModification) + ",\n")
        myRFile.write("\tmascotScore="  + str( round(-10 * math.log(1E-6 + ((float(row[8])))),2)) + ",\n")
        myRFile.write("\tproteinInformation=''"  + ",\n")
        fileName = str(row[9])
        myRFile.write("\tfileName=" + repr(fileName) + ",\n")
        myRFile.write("\trtinseconds=" + str(60 * row[6]) + "\n")

        myRFile.write("); \n \n")

    # todo run R 
    myRFile.write("save(" + myRDataName + ", file='" + myRDataName + ".RData', compress=TRUE)")


if __name__ == "__main__":
    sqlitedb=sys.argv[1]
    myquery(sqlitedb)
    print "System exit 0"

#    try:
#        t2=numpy.fromstring(row[2], dtype='float32')
#        print "size = ", sys.getsizeof(t2)
#    except:
#        pass

