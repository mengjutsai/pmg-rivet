#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os, sys
import argparse
import os.path
from os import path
import glob

def getArgs():

    parser = argparse.ArgumentParser(description='Submit condor jobs')


    parser.add_argument("--joboption", dest="joboption", action="store", default="", help="Job option to run")
    parser.add_argument("--cluster", dest="cluster", action="store", default="michigan", help="cluster to run")
    parser.add_argument("--job", dest="job", action="store", default="", help="Input for the job name")
    parser.add_argument("--applyXS", dest="applyXS", action="store_true", default=False, help="If apply Xsec")
    parser.add_argument("--script", dest="script", action="store", default="", help="script")
    parser.add_argument("--no_sub", dest="no_sub", action="store_true", default=False, help="No submission")
    parser.add_argument("--xs", dest="xs", action="store", default="589.220", help="xsec")
    parser.add_argument("--evtMax", dest="evtMax", action="store", default="500", help="evtMax")

    parser.add_argument("--input", dest="input", action="store", default="", help="input")
    parser.add_argument("--remote", dest="remote", action="store", default="", help="remote")
    parser.add_argument("--rank", dest="rank", action="store_true", default=False, help="rank for requirements")

    dict_args = vars(parser.parse_args()) # create list from args parser

    return parser.parse_args(), dict_args


def WriteCommand(dict_args):

    if args.applyXS:
        Command = 'athena  -c \"evtMax='+args.evtMax+';applyXS=True;xs='+args.xs+'\" ' + args.joboption
    else:
        Command = 'athena  -c \"evtMax='+args.evtMax+';applyXS=False\" ' + args.joboption

    return Command

def Runcluster():

    if args.cluster=="lxplus":
        path = "condor_lxplus"
    elif args.cluster=="michigan":
        path = "condor_mich"

    if args.cluster=="lxplus":
        working = "condor_lxplus_working"
    elif args.cluster=="michigan": # /atlas/data19/metsai/condor
        working = "condor_mich_working"

    return path, working


def Run(dict_args,JobFolderName="default"):

    WorkingPath=os.path.abspath(os.getcwd())

    Condor_path, Condor_working_path = Runcluster()

    output = WorkingPath+"/"+Condor_path+"/"+JobFolderName+"/"
    outputworking = WorkingPath+"/"+Condor_working_path+"/"+JobFolderName+"/"

    if path.isdir(output):
        os.system("rm -rf "+output)
    if path.isdir(outputworking):
        os.system("rm -rf "+outputworking)

    os.system("mkdir -p "+output)
    os.system("mkdir -p "+outputworking)

    os.chdir(outputworking)

    if args.input!="":
        os.system("cp "+args.input+" "+outputworking+"input.txt")

    with open(JobFolderName+'.sh', "w") as ftw:
        ftw.write("#!/bin/bash\n")
        ftw.write("path="+WorkingPath+"\n")
        ftw.write("export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase\n")
        ftw.write("source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh\n")
        ftw.write("cd $path\n")
        ftw.write("source setup.sh\n")
        ftw.write("cp library/"+args.script+" "+outputworking+"\n")
        ftw.write("cd "+outputworking+"\n")
        Command = WriteCommand(dict_args)
        ftw.write(Command+"\n")

        ftw.write("lsetup \"lcgenv -p LCG_88 $LCG_PLATFORM yoda 1.8.3\" \n")
        ftw.write("cp $path/share/yoda2root.py "+outputworking+" \n")

        outputname = args.joboption.split("/")[-1].replace(".py","")
        output_yoda = outputname+".yoda"
        output_root = outputname+".root"

        ftw.write("python yoda2root.py "+output_yoda+" \n")
        ftw.write("cp "+output_yoda+" "+output+"\n")
        ftw.write("cp "+output_root+" "+output+"\n")

        ftw.write("rm "+output_yoda+" \n")
        ftw.write("rm "+output_root+" \n")





    with open("job.sub", "w") as jobsub:

        jobsub.write("universe = vanilla\n")
        jobsub.write("executable = "+JobFolderName+'.sh'+"\n")
        jobsub.write("output = "+output+"job.out\n")
        jobsub.write("error = "+output+"job.err\n")
        jobsub.write("log = "+output+"job.log\n")
        if args.rank:
            # jobsub.write("Requirements = (Name == \"gpuatum01\") || (machine == \"gpuatum02\") || (machine == \"gpuatum03\") || (machine == \"gpuatum04\") || (machine == \"gpuatum05\") || (machine == \"pcatum02\") || (machine == \"pcatum03\") || (machine == \"pcatum04\") || (machine == \"pcatum05\") \n")
            jobsub.write("Requirements = (TARGET.Name == \"slot1@gpuatum01.dyndns.cern.ch\") \n")
        jobsub.write("queue\n")

def Submission(JobFolderName):

    os.system("pwd; ls")
    os.system('chmod +x '+JobFolderName+'.sh')
    # if args.remote!="":
    #     batch_cmd = "condor_submit -remote "+args.remote+" job.sub"
    # else:
    batch_cmd = "condor_submit job.sub"
    print ("\t > Submission of %s on the HTCaondor" % JobFolderName)
    if args.no_sub:
        pass
    else:
        os.system (batch_cmd)

def main():

    command = WriteCommand(dict_args)
    print(" >> Run command:\n%s" % command)

    JobFolderName = args.job
    Run(dict_args,JobFolderName)
    Submission(JobFolderName)

if __name__ == '__main__':
    args, dict_args = getArgs()
    main()
