// create structure with sum with your functions:

// here is an example of 8 functions with common 4 parammeters: par[1], par[4], par[5], par[6]

struct GlobalChi2PolyNfc {
   GlobalChi2PolyNfc(  ROOT::Math::IMultiGenFunction & f1,  ROOT::Math::IMultiGenFunction & f2,  ROOT::Math::IMultiGenFunction & f3,  ROOT::Math::IMultiGenFunction & f4,  ROOT::Math::IMultiGenFunction & f5,  ROOT::Math::IMultiGenFunction & f6,  ROOT::Math::IMultiGenFunction & f7,  ROOT::Math::IMultiGenFunction & f8) :
      fChi2_1(&f1), fChi2_2(&f2), fChi2_3(&f3), fChi2_4(&f4), fChi2_5(&f5), fChi2_6(&f6), fChi2_7(&f7), fChi2_8(&f8) {}

   double operator() (const double *par) const {
      double p1[10];

      p1[0] = par[0]; // Vbd 35
      p1[1] = par[1]; // a (Nc)
      p1[2] = par[2]; // Vcr
      p1[3] = par[3]; // Pgeiger
      p1[4] = par[4]; // b (Nc)
      p1[5] = par[5]; // c (Nc)
      p1[6] = par[6]; // d (Nc)
      p1[7] = par[7]; // T
      p1[8] = par[8]; // exp p0
      p1[9] = par[9]; // exp p1

      double p2[10];
      p2[0] = par[10]; // Vbd 25
      p2[1] = par[1]; // Cu * a
      p2[2] = par[11]; // Vcr
      p2[3] = par[12]; // Pgeiger
      p2[4] = par[4]; // b (Nc)
      p2[5] = par[5]; // c (Nc)
      p2[6] = par[6]; // d (Nc)
      p2[7] = par[13]; // T
      p2[8] = par[14]; // exp p0
      p2[9] = par[15]; // exp p1

      double p3[10];
      p3[0] = par[16]; // Vbd 15
      p3[1] = par[1]; // Cu * a
      p3[2] = par[17]; // Vcr
      p3[3] = par[18]; // Pgeiger
      p3[4] = par[4]; // b (Nc)
      p3[5] = par[5]; // c (Nc)
      p3[6] = par[6]; // d (Nc)
      p3[7] = par[19]; // T
      p3[8] = par[20]; // exp p0
      p3[9] = par[21]; // exp p1

      double p4[10];
      p4[0] = par[22]; // Vbd 5
      p4[1] = par[1]; // Cu * a
      p4[2] = par[23]; // Vcr
      p4[3] = par[24]; // Pgeiger
      p4[4] = par[4]; // b (Nc)
      p4[5] = par[5]; // c (Nc)
      p4[6] = par[6]; // d (Nc)
      p4[7] = par[25]; // T
      p4[8] = par[26]; // exp p0
      p4[9] = par[27]; // exp p1

      double p5[10];
      p5[0] = par[28]; // Vbd -5
      p5[1] = par[1]; // Cu * a
      p5[2] = par[29]; // Vcr
      p5[3] = par[30]; // Pgeiger
      p5[4] = par[4]; // b (Nc)
      p5[5] = par[5]; // c (Nc)
      p5[6] = par[6]; // d (Nc)
      p5[7] = par[31]; // T
      p5[8] = par[32]; // exp p0
      p5[9] = par[33]; // exp p1

      double p6[10];
      p6[0] = par[34]; // Vbd -15
      p6[1] = par[1]; // Cu * a
      p6[2] = par[35]; // Vcr
      p6[3] = par[36]; // Pgeiger
      p6[4] = par[4]; // b (Nc)
      p6[5] = par[5]; // c (Nc)
      p6[6] = par[6]; // d (Nc)
      p6[7] = par[37]; // T
      p6[8] = par[38]; // exp p0
      p6[9] = par[39]; // exp p1

      double p7[10];
      p7[0] = par[40]; // Vbd -25
      p7[1] = par[1]; // Cu * a
      p7[2] = par[41]; // Vcr
      p7[3] = par[42]; // Pgeiger
      p7[4] = par[4]; // b (Nc)
      p7[5] = par[5]; // c (Nc)
      p7[6] = par[6]; // d (Nc)
      p7[7] = par[43]; // T
      p7[8] = par[44]; // exp p0
      p7[9] = par[45]; // exp p1

      double p8[10];
      p8[0] = par[46]; // Vbd -35
      p8[1] = par[1]; // Cu * a
      p8[2] = par[47]; // Vcr
      p8[3] = par[48]; // Pgeiger
      p8[4] = par[4]; // b (Nc)
      p8[5] = par[5]; // c (Nc)
      p8[6] = par[6]; // d (Nc)
      p8[7] = par[49]; // T
      p8[8] = par[50]; // exp p0
      p8[9] = par[51]; // exp p1

      return (*fChi2_1)(p1) + (*fChi2_2)(p2) + (*fChi2_3)(p3) + (*fChi2_4)(p4) + (*fChi2_5)(p5) + (*fChi2_6)(p6) + (*fChi2_7)(p7) + (*fChi2_8)(p8);
   }

   const  ROOT::Math::IMultiGenFunction * fChi2_1;
   const  ROOT::Math::IMultiGenFunction * fChi2_2;
   const  ROOT::Math::IMultiGenFunction * fChi2_3;
   const  ROOT::Math::IMultiGenFunction * fChi2_4;
   const  ROOT::Math::IMultiGenFunction * fChi2_5;
   const  ROOT::Math::IMultiGenFunction * fChi2_6;
   const  ROOT::Math::IMultiGenFunction * fChi2_7;
   const  ROOT::Math::IMultiGenFunction * fChi2_8;

};

// create a singe fit function:
// here is a function to fit SiPM IV curve
Double_t fitIV4PolyNfc(Double_t *x, Double_t *par){

   // x[0] = Vbisa
   // Vbd = par[0];
   // Vcr = par[2];
   // T = par[7];
   //k = Pgeiger = 1. - TMath::Exp(-par[3]*(Vbias-Vbd))
   // Ipre-bd = TMath::Exp(par[8] + par[9]*Vbias);
   // Nfree-carriers = (1. + par[1]*T + par[4]*T*T)*(TMath::Exp(par[5] + par[6]*T))

   Double_t k = 1. - TMath::Exp(-par[3]*(x[0]-par[0]));
   if(k>1.0) k =1.0;
   if(k<0.0) k = 0.0;

   Double_t Paft = (x[0]-par[0])/(par[2]-par[0]);
   if(Paft>0.9999) Paft = 0.99999;
   if(Paft<0.0) Paft = 0.0;

   Double_t fitval = TMath::Exp(par[8] + par[9]*x[0]) + (1. + par[1]*par[7] + par[4]*par[7]*par[7])*(TMath::Exp(par[5] + par[6]*par[7]))*(1. + Paft/(1.-Paft) )*k*(x[0]-par[0]);

   return fitval;
}


////////////////////////////////////////////////////
// in your main program you can use a fit like this:

// define 8 function which you will use for your fit:
   TF1 *f_35 = new TF1("f_35", fitIV4PolyNfc, FitRange35C[0], FitRange35C[1], 10);           //#1
   TF1 *f_25 = new TF1("f_25", fitIV4PolyNfc, FitRange25C[0], FitRange25C[1], 10);           //#2
   TF1 *f_15 = new TF1("f_15", fitIV4PolyNfc, FitRange15C[0], FitRange15C[1], 10);           //#3
   TF1 *f_5  = new TF1("f_5" , fitIV4PolyNfc, FitRange5C[0] , FitRange5C[1] , 10);           //#4
   TF1 *f_m5 = new TF1("f_m5", fitIV4PolyNfc, FitRangem5C[0], FitRangem5C[1], 10);           //#5
   TF1 *f_m15 = new TF1("f_m5", fitIV4PolyNfc, FitRangem15C[0], FitRangem15C[1], 10);        //#6
   TF1 *f_m25 = new TF1("f_m25", fitIV4PolyNfc, FitRangem25C[0] - 0.5, FitRangem25C[1], 10); //#7
   TF1 *f_m35 = new TF1("f_m35", fitIV4PolyNfc, FitRangem35C[0]-1., FitRangem35C[1], 10);    //#8
   
   
   ROOT::Math::WrappedMultiTF1 wf35(*f_35,1);
   ROOT::Math::WrappedMultiTF1 wf25(*f_25,1);
   ROOT::Math::WrappedMultiTF1 wf15(*f_15,1);
   ROOT::Math::WrappedMultiTF1 wf5 (*f_5, 1);
   ROOT::Math::WrappedMultiTF1 wfm5(*f_m5,1);
   ROOT::Math::WrappedMultiTF1 wfm15(*f_m15,1);
   ROOT::Math::WrappedMultiTF1 wfm25(*f_m25,1);
   ROOT::Math::WrappedMultiTF1 wfm35(*f_m35,1);
   
   
   ROOT::Fit::DataOptions opt;

   // set your fit range (range35, range25, range15, range5, rangem5, rangem15, rangem25, rangem35) and data in TGRaph format (gr35C, gr25C, gr15C, gr5C, grm5C, grm15C, grm25C, grm35C) 
   ROOT::Fit::DataRange range35;
   range35.SetRange(FitRange35C[0], FitRange35C[1]);
   ROOT::Fit::BinData data35(opt,range35);
   ROOT::Fit::FillData(data35, gr35C);

   ROOT::Fit::DataRange range25;
   range25.SetRange(FitRange25C[0], FitRange25C[1]);
   ROOT::Fit::BinData data25(opt,range25);
   ROOT::Fit::FillData(data25, gr25C);

   ROOT::Fit::DataRange range15;
   range15.SetRange(FitRange15C[0], FitRange15C[1]);
   ROOT::Fit::BinData data15(opt,range15);
   ROOT::Fit::FillData(data15, gr15C);

   ROOT::Fit::DataRange range5;
   range5.SetRange(FitRange5C[0], FitRange5C[1]);
   ROOT::Fit::BinData data5(opt,range5);
   ROOT::Fit::FillData(data5, gr5C);

   ROOT::Fit::DataRange rangem5;
   rangem5.SetRange(FitRangem5C[0], FitRangem5C[1]);
   ROOT::Fit::BinData datam5(opt,rangem5);
   ROOT::Fit::FillData(datam5, grm5C);

   ROOT::Fit::DataRange rangem15;
   rangem15.SetRange(FitRangem15C[0], FitRangem15C[1]);
   ROOT::Fit::BinData datam15(opt,rangem15);
   ROOT::Fit::FillData(datam15, grm15C);

   ROOT::Fit::DataRange rangem25;
   rangem25.SetRange(FitRangem25C[0], FitRangem25C[1]);
   ROOT::Fit::BinData datam25(opt,rangem25);
   ROOT::Fit::FillData(datam25, grm25C);

   ROOT::Fit::DataRange rangem35;
   rangem35.SetRange(FitRangem35C[0], FitRangem35C[1]);
   ROOT::Fit::BinData datam35(opt,rangem35);
   ROOT::Fit::FillData(datam35, grm35C);

   // create chi-square for each function;
   
   ROOT::Fit::Chi2Function chi2_35(data35, wf35);
   ROOT::Fit::Chi2Function chi2_25(data25, wf25);
   ROOT::Fit::Chi2Function chi2_15(data15, wf15);
   ROOT::Fit::Chi2Function chi2_5(data5, wf5);
   ROOT::Fit::Chi2Function chi2_m5(datam5, wfm5);
   ROOT::Fit::Chi2Function chi2_m15(datam15, wfm15);
   ROOT::Fit::Chi2Function chi2_m25(datam25, wfm25);
   ROOT::Fit::Chi2Function chi2_m35(datam35, wfm35);

   // create a "global" (total) chi-cquare
   GlobalChi2PolyNfc globalChi2(chi2_35, chi2_25, chi2_15, chi2_5, chi2_m5, chi2_m15, chi2_m25, chi2_m35);

   ROOT::Fit::Fitter fitter;
   
   
  // set Precision and Tolerance and define the Minimazer Type
  fitter.Config().MinimizerOptions().SetPrecision(1e-14);
  fitter.Config().MinimizerOptions().SetTolerance(1e-14);
  fitter.Config().MinimizerOptions().SetMinimizerType("Minuit");

   // number of input parameters
   const int Npar = 52;

  // set starting parameters values :
  // par0 - array with input parameters for "GlobalChi2StandartF", after the fit the input parameters will be replaced by the final parameters
  fitter.Config().SetParamsSettings(Npar, par0);

  // set parameters limits, if needed, as: fitter.Config().ParSettings(i).SetLimits(low, max);, where i - paramiter number, low and max are parameter nages
  fitter.Config().ParSettings(0).SetLimits(20., 90.);
  fitter.Config().ParSettings(3).SetLimits(0.01, 3.);
  fitter.Config().ParSettings(12).SetLimits(0.01, 3.);
  fitter.Config().ParSettings(18).SetLimits(0.01, 3.);
  fitter.Config().ParSettings(24).SetLimits(0.01, 3.);
  fitter.Config().ParSettings(30).SetLimits(0.01, 3.);
  fitter.Config().ParSettings(36).SetLimits(0.01, 3.);
  fitter.Config().ParSettings(42).SetLimits(0.01, 3.);
  fitter.Config().ParSettings(48).SetLimits(0.01, 3.);
  
  // do fit
  fitter.FitFCN(Npar, globalChi2, par0, data35.Size() + data25.Size() + data15.Size() + data5.Size() + datam5.Size() + datam15.Size() + datam25.Size() + datam35.Size(), true);
  result = fitter.Result();
  result.Print(std::cout);

  // replase the input parameter by fit results
  for(int i=0; i<Npar; i++)  par0[i] = result.Value(i);
  
 //////////////////// 
