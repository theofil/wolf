#include "TCut.h"

TCut sel_basic("pt1>20 && pt2>20 && (abs(eta1)<2.4 && abs(eta2)<2.4) && abs(abs(eta1)-1.5)>0.1 && abs(abs(eta2)-1.5)>0.1 DR>0.3");
TCut sel_basic_CE("pt1>20 && pt2>20 && (abs(eta1)<1.4 && abs(eta2)<1.4) && DR>0.3");
TCut sel_basic_FW("pt1>20 && pt2>20 && (abs(eta1)<1.4 || abs(eta1)>1.6) && (abs(eta2)<1.4 || abs(eta2)>1.6) &&  (abs(eta1)<2.4 && abs(eta2)<2.4) && DR>0.3");

TCut sel_M2070("mll>20 && mll<70");
TCut sel_M120("mll>120");
TCut sel_MZ("mll>71 && mll<111");

TCut sel_ee("chid1*chid2==-1 ");
TCut sel_mumu("chid1*chid2==-4 ");
TCut sel_SF = sel_ee || sel_mumu;
TCut sel_OF("chid1*chid2==-2 ");

TCut sel_ee_trig("chid1*chid2==-1 && (trigger_bit &1 )");
TCut sel_mumu_trig("chid1*chid2==-4  && (trigger_bit &2) ");
TCut sel_SF_trig = sel_ee_trig || sel_mumu_trig;
TCut sel_OF_trig("chid1*chid2==-2 && (trigger_bit &8 || trigger_bit &4)");

TCut sel_ij5("pfJetGoodID[0]!=0 && pfJetGoodNum40>=5 ");
TCut sel_ij4("pfJetGoodID[0]!=0 && pfJetGoodNum40>=4 ");
TCut sel_ij3("pfJetGoodID[0]!=0 && pfJetGoodNum40>=3 ");
TCut sel_ij1("pfJetGoodID[0]!=0 && pfJetGoodNum40>=1 ");
TCut sel_ij2("pfJetGoodID[0]!=0 && pfJetGoodNum40>=2 ");
TCut sel_ej1("pfJetGoodID[0]!=0 && pfJetGoodNum40==1 ");
TCut sel_ej2("pfJetGoodID[0]!=0 && pfJetGoodNum40==2 ");
TCut sel_ej3("pfJetGoodID[0]!=0 && pfJetGoodNum40==3 ");
TCut sel_ej4("pfJetGoodID[0]!=0 && pfJetGoodNum40==4 ");
TCut sel_ej5("pfJetGoodID[0]!=0 && pfJetGoodNum40==5 ");

TCut sel_bij1("pfJetGoodNumBtag30>=1");
TCut sel_bij2("pfJetGoodNumBtag30>=2");
TCut sel_bej0("pfJetGoodNumBtag30==0");
TCut sel_bej1("pfJetGoodNumBtag30==1");
TCut sel_bej2("pfJetGoodNumBtag30==2");

TCut sel_met100("met[4]>100");
TCut sel_met("met[4]>100");
TCut sel_met100150("met[4]>100 && met[4]<150");
TCut sel_met150("met[4]>150");
TCut sel_mll("mll>20 && mll<70");

TCut sel_basic_gen("genPt1>20 && genPt2>20 && abs(genEta1)<1.4 && abs(genEta2)<1.4 && (genId1*genId2==-11*11 || genId1*genId2==-11*13 || genId1*genId2==-13*13)");
TCut sel_ij3_gen("genNjets>=3");
TCut sel_met_gen("genMET>100");
TCut sel_mll_gen("genMll>20 && genMll<70");
TCut sel_M2070_gen("genMll>20 && genMll<70");

TCut sel_ee_gen("genId1*genId2==-11*11 ");
TCut sel_mumu_gen("chid1*chid2==-13*13 ");
TCut sel_SF_gen = sel_ee || sel_mumu;
TCut sel_OF_gen("chid1*chid2==-11*13 ");

TCut sel_filters("passed_filters");
TCut sel_triggers("passed_triggers");
TCut sel_filter_reduced("(filter_bit==0||filter_bit==64)");

TCut fake("!(genId1Sel*genId2Sel==-11*11 || genId1Sel*genId2Sel==-13*13 || genId1Sel*genId2Sel==-11*13)");
TCut diLep("(genMID1*genMID2==-24*24 || genMID1*genGMID2==-24*24 || genGMID1*genMID2==-24*24 || genGMID1*genGMID2==-24*24)");

TCut sel_blockA("runNum<201679");
TCut sel_blockB("runNum>=201679");

TCut sel_runA("runNum>190456 && runNum<193621");
TCut sel_runB("runNum>193833 && runNum<196531");
TCut sel_runC("runNum>198022 && runNum<203742");
TCut sel_runD("runNum>203777 && runNum<208686");

TCut sel_run0("runNum>190645 && runNum<194643");
TCut sel_run1("runNum>194644 && runNum<195868");
TCut sel_run2("runNum>195915 && runNum<198955");
TCut sel_run3("runNum>198969 && runNum<200075");
TCut sel_run4("runNum>200091 && runNum<201669");
TCut sel_run5("runNum>201671 && runNum<202972");
TCut sel_run6("runNum>202973 && runNum<205344");
TCut sel_run7("runNum>205515 && runNum<206512");
TCut sel_run8("runNum>206513 && runNum<207491");
TCut sel_run9("runNum>207492 && runNum<208686");



//TCut CR_all(((((((passed_filters * passed_triggers *  (  (id1==id2)(id1==0)(trigger_bit&1) + (id1==id2)(id1==1)(trigger_bit&2) + (id1!=id2)*(trigger_bit&4|trigger_bit&8)   ))  ||!is_data))&&((abs(eta1)<1.4 && abs(eta2)<1.4 && pt1>20 && pt2>20) && abs(abs(eta1)-1.5)>0.1&&abs(abs(eta2)-1.5)>0.1 && l1l2dR>0.3))&&(((pfJetGoodNum40>=2&&pfJetGoodID[0]!=0)&&(pfJetGoodNum40>=2&&pfJetGoodID[1]!=0))&&(mll>2)))&&(mll>20))&&(pfJetGoodNum40==2&&met[4]>100&&met[4]<150))&&((id1==id2)&&(ch1*ch2<0))

