cd rivet_plugins/
rivet-build -j 24 Rivet_tttt_parton.so tttt_parton.cc
cp Rivet_tttt_parton.so ../rivet_ana/
#cp Rivet_tttt_parton.so ../library/
cd ..

# rivet-build -j 24 Rivet_parton_ttW_ttH-fixedXsec.so ttW_ttH_parton-fixedXsec.cc
# rivet-build -j 24 Rivet_parton_ttW_ttH_sh2210.so ttW_ttH_parton_sh2210.cc
# cp Rivet_parton_ttW_ttH.so Rivet_parton_ttW_ttH-fixedXsec.so Rivet_parton_ttW_ttH_sh2210.so ../rivet_ana/
# cp Rivet_parton_ttW_ttH.so Rivet_parton_ttW_ttH-fixedXsec.so Rivet_parton_ttW_ttH_sh2210.so ../library/
