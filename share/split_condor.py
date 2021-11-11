#!/usr/bin/env python

from glob import glob
import sys,os
import argparse


def getArgs():

    parser = argparse.ArgumentParser(description='Submit condor jobs')
    parser.add_argument("--input", dest="input", action="store", default="mc16_13TeV.700000.Sh_228_ttW.merge.EVNT.e7793_e7400", help="Inputfolder")
    parser.add_argument("--inputfolder", dest="inputfolder", action="store", default="/net/s3_data_home/metsai/hbsm4top/hbsm4tops-simulation/condor_mich/", help="Inputfolder")

    parser.add_argument("--nFilePerJob", dest="nFilePerJob", metavar='N', type=int, default=20, help="Number of file per job")
    parser.add_argument("--output", dest="output", action="store", default="default", help="output folder")

    # dict_args = vars(parser.parse_args()) # create list from args parser
    return parser.parse_args()

def main():

    WorkingPath=os.path.abspath(os.getcwd())

    SampleList = glob(args.inputfolder+args.input+"*/*EVNT.root")
    SampleList = sorted(SampleList)
    nFilePerJob = args.nFilePerJob

    output = WorkingPath+"/input/"+args.output+"/"

    if os.path.isdir(output):
        os.system("rm -rf "+output)
    os.system("mkdir -p "+output)

    chunks = [SampleList[x:x+nFilePerJob] for x in xrange(0, len(SampleList), nFilePerJob)]
    length = 0

    for i, sublist in enumerate(chunks):
        length = length + len(sublist)
        with open(output+'input'+str(i)+'.txt','w') as f:
            for item in sublist:
                f.write("%s\n" % item)

    print("orginial =", len(SampleList))
    print("sum of split =", length)

if __name__ == '__main__':
    args = getArgs()
    main()
