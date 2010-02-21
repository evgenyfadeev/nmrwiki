/*  CaHD_cpmg_GLY_dfh_600_v2.c
   This pulse sequence will allow t2 to perform the following experiment:
   2D C13,1H correlation map to observe exchange with variable CPMG spacing
      exchange is quantified in the Ca positions of GLY where labeled with
      50%D, i.e., we quantify the dispersion at the CHD position, while suppressing
      CH2.
   Uses two channels:
         1) 1H              - 4.73 ppm
         2) 13C             - 43.9 ppm  (so we can suppress the Cb from other residues).
         3) 15C             - 119 ppm
   Set dm = 'nnny', dmm = 'cccp' [13C decoupling during acquisition].
   Set dm2 = 'nnnn', dmm2 = 'cccp'.
   Must set phase = 1,2 for States-TPPI acquisition in t1 [13C].
   in this sequence the deuterium decoupling is performed as a series of
   pi pulses (not covering the carbon pulse).
   Proton Pulses:
   --------------
   Sit on the HDO signal at 4.774ppm
      pw      @ tpwr   : Hard pulse - max power
      pwsl    @ tpwrsl : SpinLock power during carbon CPMG.;
   Carbon Pulses:
   --------------
   Sit at 43.9 ppm with the carrier.
      pwc     @ dhpwr   : high power 13C 90 pulse at dhpwr
      pwc_cp @ dpwr_cp : pulse width of the 90s immediately flanking the CPMG period.
      pwc_ca @ d_ca     : Selective (Ca && !CO) pulses. Phase can be adjusted with phase_ca.
                          Use the RF_CO_Ca_burp_600_v3 program
      pwc_reb@ d_reb    : A selective Re-burp pulse (Ca(gly) && !Cb(!gly))
                          shape is shp_reb. Maximum length is 2.750 ms.
      1./dmf @ dpwr     : Decoupling during Acquisition
      1./dmf_co @ dpwrcodec : Decoupling during indirect acquisition
                               Pattern 'codecseq' should be a low power adiabatic decoupling,
                               e.g. a WURST-2 sequence.
      ncyc can be either even or odd :)
   Deuterium Pulses:
   -----------------
   Sit in the middle of the Glycine region (4.0)ppm
      pwd: Low power decoupling - around 350us
      pwd1: High power deuterium pulses.
   Delays
   ------
      taua: set to slightly less than 1/4J
      taub: set to exactly 1/4J
   Written D. F. Hansen on August 12 2008.
   Modified on Aug 15 2008 by DF Hansen,
     Only decouple 2H between the carbon pulses
   Modified on Aug 16th 2008 by DF Hansen,
     Select a decoupling pattern for 2H that is optimum for the tauCPMG
     delay.
     One must have 'CW.DEC','Composite.DEC','WALTZ-1.DEC', and 'WALTZ-4.DEC'
   Modified on Aug 21 2008 by DF Hansen,
     changed phase cycle (t10).
   Modified on August 29 2008 by DF Hansen
      Replace the HetCP transfer from H->C with a double refocus INEPT
      Changed phase cycle.
*/
#include <standard.h>
static int  phi1[2] = {1,3},
            phi2[4] = {1,1,3,3},
            phi3[8] = {1,1,1,1,3,3,3,3},
            phi4[4] = {0,0,0,0},
            phi6[16]= {0,0,0,0,0,0,0,0,2,2,2,2,2,2,2,2}, /* 180x in middle of CPMG */
            phi8[1] = {0},
            phi10[8]= {0,0,0,0,2,2,2,2},
            phi11[2]= {1,1},
            phi12[16]={1,1,1,1,1,1,1,1,3,3,3,3,3,3,3,3},
            rec[8]     ={0,2,2,0,0,2,2,0};
static double d2_init=0.0;
pulsesequence() {
/* DECLARE VARIABLES */
 char       fsat[MAXSTR],
            fscuba[MAXSTR],
            f1180[MAXSTR],      /* Flag to start t1 @ halfdwell             */
            shp_ca[MAXSTR],     /* CPMG pulse */
            shp_reb[MAXSTR],
            codecseq[MAXSTR],
            opt_water[MAXSTR],
            dec3pat[MAXSTR];
 int         phase,
             t1_counter,     /* used for states tppi in t1           */
             ncyc_dipfinal,
             rec_phase;
 double      tau1,           /*  t1 delay */
             taua,           /*  set to approx. 1/4JCH */
             taub,           /*  set to exactly 1/4JCH */
             pwc,            /* 90 c pulse at dhpwr             */
             tsatpwr,        /* low level 1H trans.power for presat */
             dhpwr,          /* power level for high power 13C pulses on dec1 */
              sw1,           /* sweep width in f1                    */
             dpwr_cp,        /* power level for 13C CPMG    */
             pwc_cp,         /* 13C pw for CPMG     */
             ncyc,           /* number of CPMG cycles */
             ncycmax,
              time_T2,       /* total time for CPMG trains  */
             tauCPMG,
             dhpwr2,
             pwn,
             p10,p11,p12,p13,p14,p15,p16,p17,p18,
             pw_dip,
             p_d,
             dpwr_dip,
             t_dip,
              /* Deuterium Decoupling during CPMG train */
             dpwr3_D,        /* power level for deuterium decoupling */
             dpwr3_f,        /* Fine power adjustment of the deuterium decoupling */
             pwd,            /* pw90 D for decoupling */
             pwd1,           /* Hi power D 90 */
             /* 180 pulse in the CPMG */
             pwc_ca,
             d_ca,
             phase_ca,
             /* Re-burp on Glycine */
             pwc_reb,
             d_reb,
              /* Carbonyl decoupling during indirect acquisition */
              dmf_co,
              dresco,
              dpwrcodec,
              /* Decoupling on protons during CPMG */
              pwsl,
              pwsl_f,
              pwsl_f_fst,
              tpwrsl,
              tpwrsl_f_fst=0.,
              tpwrsl_f=0.,
              numpulse, numpulsefst,
              numD,D_delay,numD_f,D_delay_f,
              time_equil,
              gt0,
              gt1,
              gt2,
              gt3,
              gt4,gt5,gt6,gt7,
              gzlvl0,
              gzlvl1,
              gzlvl2,
              gzlvl3,
              gzlvl4,gzlvl5,gzlvl6,gzlvl7,
              gzlvl10;
/* LOAD VARIABLES */
  getstr("opt_water",opt_water);
  getstr("fsat",fsat);
  getstr("f1180",f1180);
  getstr("fscuba",fscuba);
  getstr("shp_ca",shp_ca);
  getstr("codecseq",codecseq);
  getstr("shp_reb",shp_reb);
  taua    = getval("taua");
  taub    = getval("taub");
  pwc = getval("pwc");
  tpwr = getval("tpwr");
  tsatpwr = getval("tsatpwr");
  dhpwr = getval("dhpwr");
  dpwr = getval("dpwr");
  phase = (int) ( getval("phase") + 0.5);
  sw1 = getval("sw1");
  dpwr_cp  = getval("dpwr_cp");
  pwc_cp   = getval("pwc_cp");
  ncyc     = getval("ncyc");
  time_T2  = getval("time_T2");
  pwc_reb = getval("pwc_reb");
  d_reb    = getval("d_reb");
  dhpwr2 = getval("dhpwr2");
  pwn = getval("pwn");
  ncyc_dipfinal = getval("ncyc_dipfinal");
  dpwr_dip = getval("dpwr_dip");
  t_dip = getval("t_dip");
  pw_dip = getval("pw_dip");
  rec_phase= (int)(getval("rec_phase")+0.5);
  dpwr3_D = getval("dpwr3_D");
  dpwr3_f = 0.;
  pwd = getval("pwd");
  pwd1 = getval("pwd1");
  pwc_ca = getval("pwc_ca");
  d_ca = getval("d_ca");
  phase_ca = getval("phase_ca");
  dmf_co = getval("dmf_co");
  dresco = getval("dresco");
  dpwrcodec = getval("dpwrcodec");
  /* SpinLock decoupling during the CPMG train */
  pwsl    = getval("pwsl");
  tpwrsl = getval("tpwrsl");
  /* Equilibrium time before and after CPMG */
  time_equil=getval("time_equil");
  p_d = (50.0/90.0)*pw_dip;
  p10  = 4.9*p_d;
  p11  = 7.9*p_d;
  p12  = 5.0*p_d;
  p13  = 5.5*p_d;
  p14  = 0.6*p_d;
  p15  = 4.6*p_d;
  p16  = 7.2*p_d;
  p17  = 4.9*p_d;
  p18  = 7.4*p_d;
  gt0  = getval("gt0");
  gt1  = getval("gt1");
  gt2  = getval("gt2");
  gt3  = getval("gt3");
  gt4  = getval("gt4");
  gt5  = getval("gt5");
  gt6  = getval("gt6");
  gt7  = getval("gt7");
  gzlvl0 =  getval("gzlvl0");
  gzlvl1 =  getval("gzlvl1");
  gzlvl2 =  getval("gzlvl2");
  gzlvl3 =  getval("gzlvl3");
  gzlvl4 =  getval("gzlvl4");
  gzlvl5 =  getval("gzlvl5");
  gzlvl6 =  getval("gzlvl6");
  gzlvl7 =  getval("gzlvl7");
  gzlvl10=  getval("gzlvl10");
/* LOAD PHASE TABLE */
  settable(t1,2,phi1);
  settable(t2,4,phi2);
  settable(t3,8,phi3);
  settable(t4,4,phi4);
  settable(t6,16,phi6);
  settable(t7,8,rec);
  settable(t8,1,phi8);
  settable(t10,8,phi10);
  settable(t11,2,phi11);
  settable(t12,16,phi12);
/* CHECK VALIDITY OF PARAMETER RANGES */
  if ( sfrq > 750. ){
      printf(" You should change the probe. Don't use this sequence at 800MHz\n");
      abort(1);
  };
  if((dm[A] == 'y' || dm[B] == 'y' || dm[C] == 'y' )) {
       printf("incorrect dec1 decoupler flags! ");
       abort(1);
};
if((dm2[A] == 'y' || dm2[B] == 'y' || dm2[C] == 'y' || dm2[D] == 'y')) {
      printf("incorrect dec2 decoupler flags! ");
      abort(1);
};
if( tsatpwr > 6 ) {
      printf("TSATPWR too large !!! ");
      abort(1);
};
if( dpwr > 48 )      {
      printf("don't fry the probe, DPWR too large! ");
      abort(1);
};
if(dpwr_cp > 61)       {
      printf("don't fry the probe, dpwr_cp too large: < 61! ");
      abort(1);
};
if( dpwr2 > -16 )        {
      printf("don't fry the probe, DPWR2 too large! ");
      abort(1);
};
if( pw > 20.0e-6 )         {
      printf("dont fry the probe, pw too high ! ");
      abort(1);
};
if(ncyc > 40) {
      printf("ncyc_cp is too large; must be less than 41\n");
      abort(1);
};
if(time_T2 > 0.040) {
      printf("time_T2 is too large; must be less than 40 ms\n");
      abort(1);
};
ncycmax = time_T2/1e-3;
if ( ncyc > ncycmax ) {
      printf(" Maximum of ncyc is exceeded \n");
      printf(" ncyc < time_T2*1000 = %.0f \n",ncycmax);
      abort(1);
};
if(ncyc > 0) {
      tauCPMG = time_T2/(4.0*ncyc) - pwc_ca/2.;
      if(ix==1) {
         printf("nuCPMG for current experiment is (Hz): %5.3f\n",1/(4.0*(tauCPMG+pwc_ca/2.)));
      };
 } else {
      tauCPMG = time_T2/(4.0) - pwc_ca/2.;
      if(ix==1) {
         printf("nuCPMG for current experiment is (Hz): not applicable\n");
      };
 };
 if(tauCPMG - (2.0/PI)*pwc_cp - WFG_START_DELAY < 5e-6) {
      printf("tauCPMG is too small; decrease ncyc_cp\n");
      abort(1);
   };
   if(tpwrsl > 57) {
      printf("tpwrsl < 58\n");
      abort(1);
   };
   if(dpwr_dip > 56) {
      printf("dpwr_dip < 57\n");
      abort(1);
   };
   if(t_dip > 55) {
      printf("t_dip is too high; < 56\n");
      abort(1);
   };
   if(pwc_ca > 200e-6 && d_ca > 55 && shp_ca[A]=='h'){
      printf(" Too much power for CPMG-symetric inversion pulse: (pwc_ca,d_ca)\n");
      abort(1);
   };
   if(d_reb > 55 ){
      printf(" d_reb is too large; d_reb < 55dB \n");
       abort(1);
    };
    if ( pwc_reb > 2.755e-3 ){
        printf(" pwc_reb is too long; pwc_reb<2.75ms\n");
        abort(1);
    };
    if(pwc_ca > 2e-3 ) {
       printf(" pwc_ca is too long! \n");
       abort(1);
    };
    if( d_ca > 61 ){
       printf(" d_ca too high! \n");
       abort(1);
    };
    if(dpwr3_D > 51) {
       printf("Deuterium decoupling power is too high; <= 50\n");
       abort(1);
    }
    if(dpwr3 > 59) {
       printf("Deuterium pulse power is too high; <= 59\n");
       abort(1);
    }
    if(pwd>500e-6){
       printf(" Please calibrate the deuterium decoupling to pwd slightly less than 500us\n");
       abort(1);
    };
    if(gzlvl10>400){
       printf("gzlvl10 too large!\n");
       abort(1);
    };
    /* Calculate the CW proton decoupling (during CPMG) */
    numpulsefst=floor((tauCPMG - (2/PI)*pwc_cp - 2.0e-6 - POWER_DELAY)/(2.0*pwsl));
    numpulse=floor((tauCPMG+pwc_ca/2.)/(2.0*pwsl));
    /* Calculate the corresponding 90 pulse widths */
    pwsl_f_fst=(tauCPMG - (2/PI)*pwc_cp - 2.0e-6 - POWER_DELAY)/(2.0*numpulsefst);
    pwsl_f=(tauCPMG+pwc_ca/2.)/(2.0*numpulse);
    /* Calculate the corresponding powers */
    tpwrsl_f_fst=4095.0*pwsl/pwsl_f_fst;
    tpwrsl_f=4095.0*pwsl/pwsl_f;
    if(ix == 1){
       printf("numpulsesfst %5.0f, pwsl_f_fst %4.2fus, tpwrsl_f_fst
%8.2f\n",numpulsefst,pwsl_f_fst*1e6,tpwrsl_f_fst);
       printf("numpulses     %5.0f, pwsl_f    %4.2fus, tpwrsl_f
%8.2f\n",numpulse,pwsl_f*1e6,tpwrsl_f);
    };
    /* Calculate the CW Deuterium decoupling during CPMG */
    /* see how many pi pulses we can fit into tauCPMG */
    numD =floor(tauCPMG/(2.*pwd));
    numD_f=floor( (tauCPMG-2.*pwc_cp/PI-2.e-6-POWER_DELAY)/(2*pwd));
    D_delay =(tauCPMG-2.*numD*pwd)/2.;
    D_delay_f=(tauCPMG-2.*pwc_cp/PI-2.e-6-POWER_DELAY - 2.*numD_f*pwd )/2.;
    if(ncyc<0){
       numD=floor((time_T2/2.-2.*pwc_cp/PI)/(2.*pwd));
       D_delay=(time_T2/2.-2*pwc_cp/PI-2.*numD*pwd)/2.;
    };
    /* Now select the best decoupling pattern we can fit in */
    if(numD<2){
       strcpy(dec3pat,"CW");
    };
    if ( numD>1 ){ /* 90x-180y-90x */
       strcpy(dec3pat,"Composite");
    };
    if ( numD>2 ){ /* 90x-180-x 270x */
       strcpy(dec3pat,"WALTZ-1");
    };
    if ( numD>12){
        strcpy(dec3pat,"WALTZ-4");
     };
     /* Phase incrementation for hypercomplex 2D data */
     if(phase == 2) {
        tsadd(t1,1,4);
        tsadd(t11,1,4);
     };
     /* Set up f1180 tau1 = t1                   */
     tau1 = d2;
     if(f1180[A] == 'y') {
        tau1 += ( 1.0 / (2.0*sw1)) - 2.0*pwn - 2.0e-6;
        if(tau1 < 0.4e-6) {
            tau1 = 0.4e-6;
        };
     } else {
        if(tau1>2*pwn){
           tau1=tau1-2*pwn;
        }
     };
     tau1 = tau1/2.0;
     /* Calculate modifications to phases for States-TPPI acquisition        */
     if (ix == 1) {
        d2_init = d2 ;
     };
     t1_counter = (int) ( (d2-d2_init)*sw1 + 0.5 );
     /* TPPI moves artifacts to the edge of the spectrum */
     if(t1_counter % 2) {
        tsadd(t8,2,4);
        tsadd(t7,2,4);
     };
/* BEGIN ACTUAL PULSE SEQUENCE */
status(A);
  /* if ncyc == 0 apply 1H saturation at tpwrsl, tpwrsl_f for time_T2 */
  obspower(tpwrsl);
  rgpulse(pwsl,one,0.,0.);
  txphase(zero);
  obspwrf(tpwrsl_f);
  obsunblank();
  xmtron();
  if(ncyc==0) {
     delay(time_T2);
  };
  if( (ni-1)/sw1-d2>0.){
      delay((ni-1)/sw1-d2);
  };
  xmtroff();
  obsblank();
  obspwrf(4095.0);
  rgpulse(pwsl,three,0.,0.);
  obspower(tpwrsl);
  txphase(zero);
  delay(2.0e-6);
  rgradient('z',gzlvl7);
  delay(gt7);
  rgradient('z',0.0);
  delay(200.0e-6);
  obspower(tsatpwr);        /* Set transmitter power for 1H presaturation */
  decpower(dhpwr);          /* Set Dec1 power for 13C high power pulses */
  dec2power(dhpwr2);         /* Set Dec2 power for 15N hard pulses */
  dec3power(dpwr3);         /* Set Dec3 power for 2D decoupling   */
/* Presaturation Period */
status(B);
  if (fsat[0] == 'y')     {
     delay(2.0e-5);
     rgpulse(d1,zero,2.0e-6,2.0e-6); /* presaturation */
     obspower(tpwr);     /* Set transmitter power for hard 1H pulses */
     delay(2.0e-5);
     if(fscuba[0] == 'y') {
        delay(2.2e-2);
        rgpulse(pw,zero,2.0e-6,0.0);
        rgpulse(2*pw,one,2.0e-6,0.0);
        rgpulse(pw,zero,2.0e-6,0.0);
        delay(2.2e-2);
     };
  } else {
     delay(d1);
  };
  obspower(tpwr);            /* Set transmitter power for hard 1H pulses */
  txphase(zero);
  dec2phase(zero);
  decphase(zero);
  delay(1.0e-5);
/* BEGIN ACTUAL SEQUENCE */
status(C);
  /* Let the purge element phase-cycle together with initial DIPSI
      proton phase */
  /* hlv(ct,v3); */
  hlv(ct,v2);
  mod2(v2,v1);
  hlv(ct,v2);
  hlv(v2,v3);
  mod2(v3,v2);
  /* Turn receiver off    and hold the lock */
  rcvroff();
  lk_hold();
  decpower(d_reb);
  delay(20.0e-6);
  rgpulse(pw,zero,0.0,0.0);
  delay(taub
             -WFG2_START_DELAY - 202e-6 -2*GRADIENT_DELAY - gt1 -pwc_reb/2. );
  delay(2.0e-6);
  rgradient('z',gzlvl1);
  delay(gt1);
  rgradient('z',0.0);
  delay(200.0e-6);
  simshaped_pulse("hard",shp_reb,2.*pw,pwc_reb,one,zero,0.,0.);
  txphase(t3);
  delay(taub - pwc_reb/2. - WFG2_STOP_DELAY - 202.e-6 -2*GRADIENT_DELAY - gt1 );
  delay(2.0e-6);
  rgradient('z',gzlvl1);
  delay(gt1);
  rgradient('z',0.0);
  delay(200.0e-6);
  rgpulse(pw,t3,0.,0.);
  decpower(dhpwr);
delay(2.0e-6);
rgradient('z',gzlvl2);
delay(gt2);
rgradient('z',0.0);
delay(200.0e-6);
decrgpulse(pwc,t2,0.,0.);
decpower(d_reb);
delay(taub - POWER_DELAY
          -WFG2_START_DELAY - 202e-6 -2*GRADIENT_DELAY - gt3 -pwc_reb/2. );
delay(2.0e-6);
rgradient('z',gzlvl3);
delay(gt3);
rgradient('z',0.0);
delay(200.0e-6);
simshaped_pulse("hard",shp_reb,2.*pw,pwc_reb,one,zero,0.,0.);
decpower(dhpwr);
delay(taub - POWER_DELAY
            -pwc_reb/2. - WFG2_STOP_DELAY - 202.e-6 -2*GRADIENT_DELAY - gt3 -2.*POWER_DELAY );
delay(2.0e-6);
rgradient('z',gzlvl3);
delay(gt3);
rgradient('z',0.0);
delay(200.0e-6);
decrgpulse(pwc,zero,0.,0.);
decpower(dhpwr);
obspower(tpwrsl-3.);
txphase(one);
delay(2.0e-6);
rgradient('z',gzlvl4);
delay(gt4);
rgradient('z',0.0);
delay(200.0e-6);
dec3power(dpwr3);
dec3rgpulse(pwd1,one,4.0e-6,0.);
dec3power(dpwr3_D);
dec3phase(zero);
txphase(one);
xmtron();
delay(time_equil/2.);
xmtroff();
txphase(one);
xmtron();
delay(time_equil/2.);
xmtroff();
obspower(tpwrsl); /* lower power */
obspwrf(tpwrsl_f_fst); /* Set fine power for the last tcpmg delay */
rlpower(dpwr_cp,DODEV);
decrgpulse(pwc_cp,t10,0.,0.);
decphase(one);
decpower(d_ca);
dec3power(dpwr3_D);
/* CPMG BEGIN */
if(ncyc<0){
  txphase(one);
  xmtron();
   delay(D_delay);
   dec3phase(zero);
   dec3unblank();
   dec3prgon(dec3pat,pwd,90.);
   dec3on();
     delay(time_T2/2.-2*pwc_cp/PI-2*D_delay);
   dec3off();
   dec3prgoff();
   dec3blank();
   delay(D_delay);
};
if (ncyc > 0)      {
   txphase(one);
   xmtron();
   dec3phase(zero);
   delay(D_delay_f);
   dec3unblank();
   dec3prgon(dec3pat,pwd,90.);
   dec3on();
   delay(tauCPMG - (2/PI)*pwc_cp - 2.0e-6-POWER_DELAY - 2.*D_delay_f );
   dec3off();
   dec3prgoff();
   dec3blank();
   delay(D_delay_f-WFG_START_DELAY);
   xmtroff();
   obspwrf(tpwrsl_f);
   xmtron();
   decshaped_pulse(shp_ca,pwc_ca,one,0.,0.);
   delay(D_delay-WFG_STOP_DELAY);
   dec3unblank();
   dec3prgon(dec3pat,pwd,90.);
   dec3on();
     delay(tauCPMG-2.*D_delay);
   dec3off();
   dec3prgoff();
   dec3blank();
   delay(D_delay);
};
if(ncyc > 1)     {
   initval((ncyc-1),v4);
   loop(v4,v5);
   delay(D_delay);
   dec3unblank();
   dec3prgon(dec3pat,pwd,90.);
   dec3on();
     delay(tauCPMG-2.*D_delay);
   dec3off();
   dec3prgoff();
   dec3blank();
   delay(D_delay-WFG_START_DELAY);
   decshaped_pulse(shp_ca,pwc_ca,one,0.,0.);
   delay(D_delay-WFG_STOP_DELAY);
   dec3unblank();
   dec3prgon(dec3pat,pwd,90.);
   dec3on();
     delay(tauCPMG-2.*D_delay);
   dec3off();
   dec3prgoff();
   dec3blank();
   delay(D_delay);
   endloop(v5);
};
/* CPMG END */
/* Inversion pulse ( off-resonance effects etc. ) */
/* BEGIN */
if (ncyc!=0){
   xmtroff();
   rlpower(d_ca,DODEV);
   decphase(t6);
   initval(1.0,v3);
   decstepsize(phase_ca);
   dcplrphase(v3);
   delay(WFG_STOP_DELAY);
   decshaped_pulse(shp_ca,pwc_ca,t6,0.,0.);
   delay(WFG_START_DELAY);
   dcplrphase(zero);
   rlpower(dpwr_cp,DODEV);
   decphase(one);
   txphase(one);
   delay(0.5e-6);
   xmtron();
} else {
   rlpower(d_ca,DODEV);
   decphase(t6);
   initval(1.0,v3);
   decstepsize(phase_ca);
   dcplrphase(v3);
   delay(WFG_STOP_DELAY+0.5e-6);
   decshaped_pulse(shp_ca,pwc_ca,t6,0.,0.);
   delay(WFG_START_DELAY);
   dcplrphase(zero);
   rlpower(dpwr_cp,DODEV);
   decphase(one);
   txphase(one);
   delay(0.5e-6);
};
/* Inversion pulse END */
/* CPMG - BEGIN */
if(ncyc > 1) {
   initval((ncyc-1),v13);
   loop(v13,v14);
   delay(D_delay);
   dec3unblank();
   dec3prgon(dec3pat,pwd,90.);
   dec3on();
     delay(tauCPMG-2.*D_delay);
   dec3off();
   dec3prgoff();
   dec3blank();
   delay(D_delay-WFG_START_DELAY);
     decshaped_pulse(shp_ca,pwc_ca,one,0.,0.);
   delay(D_delay-WFG_STOP_DELAY);
   dec3unblank();
   dec3prgon(dec3pat,pwd,90.);
   dec3on();
     delay(tauCPMG-2.*D_delay);
   dec3off();
   dec3prgoff();
   dec3blank();
   delay(D_delay);
   endloop(v14);
};
if(ncyc > 0) {
   delay(D_delay);
   dec3unblank();
   dec3prgon(dec3pat,pwd,90.);
   dec3on();
     delay(tauCPMG-2.*D_delay);
   dec3off();
   dec3prgoff();
   dec3blank();
   delay(D_delay-WFG_START_DELAY);
   decshaped_pulse(shp_ca,pwc_ca,one,0.,0.);
   xmtroff();
   obspwrf(tpwrsl_f_fst);
   txphase(one);
   xmtron();
   delay(D_delay_f-WFG_STOP_DELAY);
   dec3unblank();
   dec3prgon(dec3pat,pwd,90.);
   dec3on();
     delay(tauCPMG - (2/PI)*pwc_cp - 2.0e-6-POWER_DELAY -2*D_delay_f );
   dec3off();
   dec3prgoff();
   dec3blank();
   delay(D_delay_f);
   xmtroff();
};
if(ncyc<0){
    txphase(one);
    xmtron();
    delay(D_delay);
    dec3unblank();
    dec3prgon(dec3pat,pwd,90.);
    dec3on();
      delay(time_T2/2.-pwc_cp*2./PI-2*D_delay);
    dec3off();
    dec3prgoff();
    dec3blank();
    delay(D_delay);
    xmtroff();
 };
decpower(dpwr_cp);
obspwrf(4095.0);
decrgpulse(pwc_cp,zero,2.0e-6,0.0);
decpower(dhpwr);
obspower(tpwrsl-3);
txphase(one);
xmtron();
delay(time_equil/2.);
xmtroff();
txphase(one);
xmtron();
delay(time_equil/2.);
xmtroff();
obspower(tpwr);
txphase(one);
dec3power(dpwr3);
dec3rgpulse(pwd1,three,2.e-6,0.);
dec3power(dpwr3_D);
dec3phase(zero);
decphase(t1);
/* purge deuterium before starting the lock
    (for probes where we do not align 2H magnetization before spin lock) */
delay(2.0e-6);
rgradient('z',gzlvl0);
delay(gt0);
rgradient('z',0.0);
delay(500.0e-6);
rgpulse(pw,zero,0.,0.);
delay(10e-6);
/* Carbon excitation pulse */
if(opt_water[A]=='y'){
   decrgpulse(pwc,t11,0.,0.);
} else {
   decrgpulse(pwc,t1,0.,0.);
};
decpower(d_ca);
decphase(zero);
delay(taub - POWER_DELAY-202e-6 -2*GRADIENT_DELAY - gt5 -pwc_ca/2. - WFG2_START_DELAY );
delay(2.0e-6);
rgradient('z',gzlvl5);
delay(gt5);
rgradient('z',0.0);
delay(200.0e-6);
simshaped_pulse("hard",shp_ca,2*pw,pwc_ca,one,zero,0.,0.);
delay(taub -pwc_ca/2. -WFG2_STOP_DELAY
            - 202e-6 -2*GRADIENT_DELAY - gt5 - 2.*pw - 3*POWER_DELAY - pwsl
            - 2*POWER_DELAY -pwd1 - 4.0e-6 -POWER_DELAY );
                                          /*     \- from decpower after 2*tau1 */
delay(2.0e-6);
rgradient('z',gzlvl5);
delay(gt5);
rgradient('z',0.0);
delay(200.0e-6);
delay(2.*pw);
dec2power(dhpwr2);
decpower(dpwrcodec);
/* Turn on proton spin lock */
obspower(tpwrsl);
rgpulse(pwsl,t4,0.,0.);
txphase(t12);
xmtron();
/* turn on carbonyl decoupling */
decprgon(codecseq,1.0/dmf_co,dresco);
decon();
/* turn on deuterium decoupling */
dec3phase(one);
dec3power(dpwr3);
dec3rgpulse(pwd1,one,4.0e-6,0.0);
dec3phase(zero);
dec3unblank();
dec3power(dpwr3_D);
dec3prgon("waltz16",pwd,90.);
dec3on();
  delay(tau1);
  dec2rgpulse(2.*pwn,zero,0.,0.);
  delay(tau1);
  decoff();
  decprgoff();
  decpower(dhpwr);
  decrgpulse(pwc,t8,0.,0.);
  /* Turn off deuterium decoupling */
  dec3off();
  dec3prgoff();
  dec3blank();
  dec3phase(three);
  dec3power(dpwr3);
  dec3rgpulse(pwd1,three,4.0e-6,0.0);
  xmtroff();
  obspower(tpwr);
  decpower(d_reb);
  delay(taub - 2*POWER_DELAY
            -WFG2_START_DELAY - 202e-6 -2*GRADIENT_DELAY - gt6 -pwc_reb/2. );
  delay(2.0e-6);
  rgradient('z',gzlvl6);
  delay(gt6);
  rgradient('z',0.0);
  delay(200.0e-6);
  simshaped_pulse("hard",shp_reb,2.*pw,pwc_reb,one,zero,0.,0.);
  delay(taub - pwc_reb/2. - WFG2_STOP_DELAY - 202.e-6 -2*GRADIENT_DELAY - gt6 -2.*POWER_DELAY );
  delay(2.0e-6);
  rgradient('z',gzlvl6);
  delay(gt6);
  rgradient('z',0.0);
  delay(200.0e-6-4.*pw/PI);
  rgpulse(pw,zero,0.,0.);
  decpower(dpwr);    /* Set power for decoupling */
  dec2power(dpwr2); /* Set power for decoupling */
  lk_sample();       /* Turn on the lock */
  /* BEGIN ACQUISITION */
status(D);
  setreceiver(t7);
}
