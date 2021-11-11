cd rivet_plugins/
rivet-build -j 24 Rivet_tttt_event.so tttt_event.cc
cp Rivet_tttt_event.so  ../rivet_ana/
cp Rivet_tttt_event.so  ../library/
cd ..



# cd rivet_plugins/
# rivet-build -j 24 Rivet_ttW_ttH_analysis.so ttw_ttH-ml.cc
# rivet-build -j 24 Rivet_ttW_ttH_powheg_analysis.so ttw_ttH-ml_powheg.cc
# cp Rivet_ttW_ttH_analysis.so Rivet_ttW_ttH_powheg_analysis.so ../rivet_ana/
# cp Rivet_ttW_ttH_analysis.so Rivet_ttW_ttH_powheg_analysis.so ../library/
# cd ..


# rivet-build -j 24 Rivet_ttW_ttH_analysis-fixedXsec.so ttw_ttH-ml-fixedXsec.cc
# rivet-build -j 24 Rivet_ttW_ttH_analysis_sh2210.so ttw_ttH-ml_sh2210.cc
# rivet-build -j 24 Rivet_MC_ttW_inclusive.so MC_ttW_inclusive.cc
# rivet-build -j 24 Rivet_MC_ttW_2SameSignLeptons.so MC_ttW_2SameSignLeptons.cc

# cp Rivet_ttW_ttH_analysis.so Rivet_ttW_ttH_analysis-fixedXsec.so Rivet_ttW_ttH_analysis_sh2210.so ../rivet_ana/
# cp Rivet_ttW_ttH_analysis.so Rivet_ttW_ttH_analysis-fixedXsec.so Rivet_ttW_ttH_analysis_sh2210.so ../library/
# cp Rivet_MC_ttW_inclusive.so Rivet_MC_ttW_2SameSignLeptons.so  ../rivet_ana/
# cp Rivet_MC_ttW_inclusive.so Rivet_MC_ttW_2SameSignLeptons.so  ../library/
