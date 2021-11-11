# pmg-rivet

First time, initialize the setup
```
source env.sh
```

Every time login later,
```
source setup.sh
source init.sh
source init.sh
```

To rebuild the source files,
```
# Rebuild the event level
source rebuild.sh
# Rebuild the parton level
source rebuild_parton.sh
```

To run the example, see the `script/run_event_tttt.sh` and `script/run_parton_tttt.sh`. You need to change the input and source library with your framework path.

There are several steps,
1. Use `share/split_condor.py` to get the input list from a set of job folders. The list is later used as input list for the rivet running.
2. Use `share/submit.py` to submit the job to run through the EVNT files to generate the ROOT output histogram file. See example `script/run_event_tttt.sh` and `script/run_parton_tttt.sh`.
3. Once you produce several output root and yoda files. You need to merge the yoda files. See example in `script/merge.sh`.
4. Convert the merged yoda file into ROOT file. With the following command,
```
python yoda2root.py example.yoda
```
It will convert yoda file to ROOT file with same name. (The output will become `example.root`)
5. Then finally, you need to write the script to draw the histograms from ROOT file.
