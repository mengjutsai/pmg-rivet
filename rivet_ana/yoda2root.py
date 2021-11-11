from array import array
import ROOT as rt
import yoda,sys

fName = sys.argv[1]
yodaAOs = yoda.read(fName)
rtFile = rt.TFile(fName[:fName.find('.yoda')] + '.root', 'recreate')
# print(yodaAOs)
for name in yodaAOs:
  new_name = name.split("/")[-1]
  # print("new_name = %s" % new_name)

  if '[' in new_name and ']' in new_name:
      variable = new_name.split("[")[0]
      variation = new_name.split("[")[1].replace(']','')
      # print("variable = %s, variation = %s" % (variable,variation))
      stored_name = variable+"_"+variation
  else:
      stored_name = new_name

  yodaAO = yodaAOs[name];  rtAO = None
  if 'Histo1D' in str(yodaAO):
    rtAO = rt.TH1D(stored_name, '', yodaAO.numBins(), array('d', yodaAO.xEdges()))
    rtAO.Sumw2(); rtErrs = rtAO.GetSumw2()
    for i in range(rtAO.GetNbinsX()):
      rtAO.SetBinContent(i + 1, yodaAO.bin(i).sumW())
      rtErrs.AddAt(yodaAO.bin(i).sumW2(), i+1)
  elif 'Scatter2D' in str(yodaAO):
    rtAO = rt.TGraphAsymmErrors(yodaAO.numPoints())
    for i in range(yodaAO.numPoints()):
      x = yodaAO.point(i).x(); y = yodaAO.point(i).y()
      xLo, xHi = yodaAO.point(i).xErrs()
      yLo, yHi = yodaAO.point(i).yErrs()
      rtAO.SetPoint(i, x, y)
      rtAO.SetPointError(i, xLo, xHi, yLo, yHi)
  else:
    continue
  rtAO.Write(stored_name)
rtFile.Close()
