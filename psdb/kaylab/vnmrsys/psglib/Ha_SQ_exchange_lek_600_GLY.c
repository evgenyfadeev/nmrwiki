/*    Ha_SQ_exchange_lek_600_GLY.c
      This pulse sequence will allow one to perform the following experiment:
      2D C13,1H correlation map to observe exchange with variable CPMG spacing
      Used to record 1H dispersions of Ha
      Uses two channels:
           1) 1H                - 4.73 ppm or methyl region
           2) 13C               - 58 ppm (Ca region)
      Set dm = 'nnny', dmm = 'cccp' [13C decoupling during acquisition].
      Set dm2 = 'nnnn', dmm2 = 'cccp'.
      Must set phase = 1,2 for States-TPPI acquisition in t1 [13C].
      Written by L.E.Kay on Aug 6, 2008 from Ha_SQ_exchange_lek_600_v2.c
        --- sequence only for GLY residues (does not attempt to suppress CH2, since this would
             require further 1/JHC delays during which signal relaxes.
      Use a 3 ms reburp pulse (at 600 MHz) centered at 45 ppm (GLY).
      Make sure that the 2H hard pulses are turned off at the 800!
      Set MQ_flg to y if you want MQ evolution during t1
      Modified by D.F Hansen on August 8 2008
        Include a filter to remove CH2 in 1/(2*J).
*/
#include <standard.h>
static int    phi1[4]  = {0,0,2,2},
              phi2[2]  = {0,2},
              phi4[8]  = {0,0,0,0,2,2,2,2},
              rec[4]   = {0,2,2,0};
static double d2_init=0.0;
pulsesequence()
{
/* DECLARE VARIABLES */
  char        fsat[MAXSTR],
              fscuba[MAXSTR],
              f1180[MAXSTR],      /* Flag to start t1 @ halfdwell           */
              sh_reb[MAXSTR],
              codec[MAXSTR],
              MQ_flg[MAXSTR],
              filter_flg[MAXSTR];
  int          phase,
               t1_counter;    /* used for states tppi in t1          */
  double       tau1,          /* t1 delay */
               taua,          /* set to exactly 1/4JCH */
               tsatpwr,       /* low level 1H trans.power for presat */
               sw1,           /* sweep width in f1                   */
               tpwr_cp,       /* power level for 1H CPMG     */
               pw_cp,         /* 1H pw for CPMG      */
               ncyc_cp,       /* number of CPMG cycles */
               time_T2,       /* total time for CPMG trains   */
               tau_cpmg,
               dhpwr,
               pwc,
               dmf_co,
               dpwr_co,
               dresco,
              gt0,
              gt1,
              gt2,
              gt3,
              gt4,
              gt5,
              gt6,
              gzlvl0,
              gzlvl1,
              gzlvl2,
              gzlvl3,
              gzlvl4,
              gzlvl5,
              gzlvl6,
              tpwr,
              pw,
              d_reb,
              pwc_reb,
              dpwr3_D,
              pwd,
              pwd1,
              tau_eq,
              pwn,
              dhpwr2;
/* LOAD VARIABLES */
  getstr("fsat",fsat);
  getstr("f1180",f1180);
  getstr("fscuba",fscuba);
  getstr("sh_reb",sh_reb);
  getstr("codec",codec);
  getstr("MQ_flg",MQ_flg);
  getstr("filter_flg",filter_flg);
  taua    = getval("taua");
  tpwr = getval("tpwr");
  tsatpwr = getval("tsatpwr");
  dpwr = getval("dpwr");
  phase = (int) ( getval("phase") + 0.5);
  sw1 = getval("sw1");
  tpwr_cp  = getval("tpwr_cp");
  pw_cp =  getval("pw_cp");
  ncyc_cp  = getval("ncyc_cp");
  time_T2  = getval("time_T2");
  dhpwr = getval("dhpwr");
  pwc = getval("pwc");
  pwn = getval("pwn");
  dhpwr2 = getval("dhpwr2");
  dmf_co = getval("dmf_co");
  dpwr_co = getval("dpwr_co");
  dresco = getval("dresco");
  gt0  = getval("gt0");
  gt1  = getval("gt1");
  gt2  = getval("gt2");
  gt3  = getval("gt3");
  gt4  = getval("gt4");
  gt5  = getval("gt5");
  gt6 = getval("gt6");
  gzlvl0  = getval("gzlvl0");
  gzlvl1  = getval("gzlvl1");
  gzlvl2  = getval("gzlvl2");
  gzlvl3  = getval("gzlvl3");
  gzlvl4  = getval("gzlvl4");
  gzlvl5  = getval("gzlvl5");
  gzlvl6  = getval("gzlvl6");
  tpwr = getval("tpwr");
  pw = getval("pw");
  d_reb = getval("d_reb");
  pwc_reb = getval("pwc_reb");
  dpwr3_D = getval("dpwr3_D");
  pwd = getval("pwd");
  pwd1 = getval("pwd1");
  tau_eq = getval("tau_eq");
/* LOAD PHASE TABLE */
  settable(t1,4,phi1);
  settable(t2,2,phi2);
  settable(t4,8,phi4);
  settable(t5,4,rec);
/* CHECK VALIDITY OF PARAMETER RANGES */
    if((dm[A] == 'y' || dm[B] == 'y' || dm[C] == 'y' ))
    {
         printf("incorrect dec1 decoupler flags! ");
         abort(1);
    }
    if((dm2[A] == 'y' || dm2[B] == 'y' || dm2[C] == 'y' || dm2[D] == 'y'))
    {
         printf("incorrect dec2 decoupler flags! ");
         abort(1);
    }
    if( tsatpwr > 6 )
    {
         printf("TSATPWR too large !!!  ");
         abort(1);
    }
    if( dpwr > 48 )
    {
         printf("don't fry the probe, DPWR too large!  ");
         abort(1);
    }
    if(tpwr_cp > 62)
    {
         printf("don't fry the probe, tpwr_cp too large: < 62! ");
         abort(1);
    }
    if(pw_cp < 9.5e-6) {
         printf("pw_cp is too low; > 9.5us\n");
         abort(1);
    }
    if( dpwr2 > -16 )
    {
         printf("don't fry the probe, DPWR2 too large!  ");
         abort(1);
    }
      if( pw > 20.0e-6 )
      {
            printf("dont fry the probe, pw too high ! ");
            abort(1);
      }
      if(gt1 > 3e-3 || gt2 > 3e-3 || gt3 > 3e-3 || gt4 > 3e-3
             || gt5 > 3e-3 || gt6 > 3e-3)
      {
            printf("gradients on for too long. Must be < 3e-3 \n");
            abort(1);
      }
      if(ncyc_cp > 80) {
            printf("ncyc_cp is too large; must be less than 81\n");
            abort(1);
      }
      if(time_T2 > .080) {
            printf("time_T2 is too large; must be less than 80 ms\n");
            abort(1);
      }
      if(ncyc_cp > 0) {
            tau_cpmg = time_T2/(4.0*ncyc_cp) - pw_cp;
            if(ix==1)
              printf("nuCPMG for curent experiment is (Hz): %5.3f\n",1/(4.0*(tau_cpmg+pw_cp)));
        }
        else {
            tau_cpmg = time_T2/(4.0) - pw_cp;
            if(ix==1)
              printf("nuCPMG for curent experiment is (Hz): not applicable");
        }
      if(tau_cpmg + pw_cp < 125e-6) {
           printf("tau_cpmg is too small; decrease ncyc_cp\n");
           abort(1);
      }
      if(dpwr_co > 42) {
           printf("dpwr_co is too high; < 42\n");
           abort(1);
    }
      if(dpwr3_D > 51) {
           printf("dpwr3_D is too high; < 52\n");
           abort(1);
      }
      if(dpwr3 > 59) {
           printf("dpwr3 is too high; < 60\n");
           abort(1);
    }
    if(ix==1)
        printf("If at 800 turn dpwr3=-16, pwd1=0\n");
/*    Phase incrementation for hypercomplex 2D data */
      if (phase == 2)
          tsadd(t1,1,4);
/*    Set up f1180    tau1 = t1               */
      tau1 = d2;
   if(MQ_flg[A] == 'n')
      tau1 = tau1 - 4.0/PI*pwc - POWER_DELAY - PRG_START_DELAY
                - 2.0*pw - 2.0*pwn - PRG_STOP_DELAY - POWER_DELAY
                - 4.0e-6;
    else
      tau1 = tau1 - 4.0/PI*pwc - POWER_DELAY - PRG_START_DELAY
              - 2.0*pw - 2.0*pwn - PRG_STOP_DELAY - POWER_DELAY
              - 4.0e-6;
     if(f1180[A] == 'y') {
          tau1 += ( 1.0 / (2.0*sw1));
          if(tau1 < 0.4e-6) tau1 = 0.4e-6;
     }
     if(tau1 < 0.4e-6)
          tau1 = 0.4e-6;
       tau1 = tau1/2.0;
/* Calculate modifications to phases for States-TPPI acquisition              */
   if( ix == 1) d2_init = d2 ;
   t1_counter = (int) ( (d2-d2_init)*sw1 + 0.5 );
   if(t1_counter % 2) {
        tsadd(t1,2,4);
        tsadd(t5,2,4);
     }
/* BEGIN ACTUAL PULSE SEQUENCE */
status(A);
   rlpower(tsatpwr,TODEV);       /* Set transmitter power for 1H presaturation */
   rlpower(dhpwr,DODEV);         /* Set Dec1 power for 13C pulses          */
   rlpower(dhpwr2,DO2DEV);       /* Set Dec2 power for 15N pulses       */
   obsoffset(tof);
/* Presaturation Period */
status(B);
   if (fsat[0] == 'y')
   {
          delay(2.0e-5);
          rgpulse(d1,zero,2.0e-6,2.0e-6); /* presaturation */
          rlpower(tpwr,TODEV);    /* Set transmitter power for hard 1H pulses */
          delay(2.0e-5);
          if(fscuba[0] == 'y')
          {
                  delay(2.2e-2);
                  rgpulse(pw,zero,2.0e-6,0.0);
                  rgpulse(2*pw,one,2.0e-6,0.0);
                  rgpulse(pw,zero,2.0e-6,0.0);
                  delay(2.2e-2);
          }
   }
   else
   {
     delay(d1);
   }
   rlpower(tpwr,TODEV);               /* Set transmitter power for 1H CPMG pulses */
   txphase(zero);
   dec2phase(zero);
   decphase(zero);
   delay(1.0e-5);
/* Begin Pulses */
status(C);
   rcvroff();
   delay(20.0e-6);
   decrgpulse(pwc,zero,4.0e-6,0.0);
delay(2.0e-6);
rgradient('z',gzlvl1);
delay(gt1);
rgradient('z',0.0);
delay(250.0e-6);
rgpulse(pw,zero,0.0,0.0);
decpower(d_reb);
delay(2.0e-6);
rgradient('z',gzlvl2);
delay(gt2);
rgradient('z',0.0);
delay(150.0e-6);
if(filter_flg[A] == 'y')
   delay(taua - POWER_DELAY - gt2 - 152e-6
     - WFG2_START_DELAY - 0.5*pwc_reb - 4.0/PI*pw);
else
   delay(taua - POWER_DELAY - gt2 - 152e-6
     - WFG2_START_DELAY - 0.5*pwc_reb);
simshaped_pulse("hard",sh_reb,2.0*pw,pwc_reb,zero,zero,0.0,0.0);
txphase(one);
decpower(dhpwr);
decphase(t4);
delay(taua - 0.5*pwc_reb - WFG2_STOP_DELAY - POWER_DELAY - gt2 - 152e-6 );
delay(2.0e-6);
rgradient('z',gzlvl2);
delay(gt2);
rgradient('z',0.0);
delay(150.0e-6);
if(filter_flg[A] == 'n')
  rgpulse(pw,one,0.0,0.0);
if(filter_flg[A] == 'y') {
decrgpulse(pwc,t4,0.,0.);
decpower(d_reb);
decphase(zero);
delay(2.0e-6);
rgradient('z',gzlvl0);
delay(gt0);
rgradient('z',0.0);
delay(150.0e-6);
delay(taua - POWER_DELAY - gt0 - 152e-6
     - WFG2_START_DELAY - 0.5*pwc_reb);
simshaped_pulse("hard",sh_reb,2.0*pw,pwc_reb,zero,zero,0.0,0.0);
txphase(one);
decpower(dhpwr);
decphase(t4);
delay(taua - 0.5*pwc_reb - WFG2_STOP_DELAY - POWER_DELAY - gt0 - 152e-6 );
delay(2.0e-6);
rgradient('z',gzlvl0);
delay(gt0);
rgradient('z',0.0);
delay(150.0e-6);
decrgpulse(pwc,t4,0.0,0.0);
rgpulse(pw,one,0.0,0.0);
}
  decphase(t1);
  delay(2.0e-6);
  rgradient('z',gzlvl3);
  delay(gt3);
  rgradient('z',0.0);
  delay(250.0e-6);
  /* turn on 2H decoupling */
  dec3phase(one);
  dec3power(dpwr3);
  dec3rgpulse(pwd1,one,4.0e-6,0.0);
  dec3phase(zero);
  dec3unblank();
  dec3power(dpwr3_D);
  dec3prgon(dseq3,pwd,dres3);
  dec3on();
  /* turn on 2H decoupling */
  if(MQ_flg[A] == 'y') {
     rgpulse(pw,zero,2.0e-6,0.0);
     delay(2.0*pwn - PRG_START_DELAY - PRG_STOP_DELAY);
  }
  decrgpulse(pwc,t1,4.0e-6,0.0);
  decphase(zero);
  /* 13CO decoupling on */
  decpower(dpwr_co);
  decprgon(codec,1.0/dmf_co,dresco);
  decon();
  /* 13CO decoupling on */
  delay(tau1);
  rgpulse(2.0*pw,zero,0.0,0.0);
  dec2rgpulse(2.0*pwn,zero,0.0,0.0);
  delay(tau1);
  /* 13CO decoupling off */
  decoff();
  decprgoff();
  /* 13CO decoupling off */
  decpower(dhpwr);
  decrgpulse(pwc,zero,4.0e-6,0.0);
  if(MQ_flg[A] == 'y')
    rgpulse(pw,zero,0.0,0.0);
  /* turn off decoupling */
  dec3off();
  dec3prgoff();
  dec3blank();
  dec3phase(three);
  dec3power(dpwr3);
  dec3rgpulse(pwd1,three,4.0e-6,0.0);
  /* turn off decoupling */
  obspower(tpwr_cp);
  if(MQ_flg[A] == 'n') {
    delay(2.0e-6);
    rgradient('z',gzlvl4);
    delay(gt4);
    rgradient('z',0.0);
    delay(250.0e-6);
  }
  else {
    delay(2.0e-6);
    rgradient('z',-1.0*gzlvl4);
    delay(gt4);
    rgradient('z',0.0);
    delay(250.0e-6);
}
    /* now include a delay to allow the spin system to equilibrate */
    delay(tau_eq);
    rgpulse(pw_cp,t2,4.0e-6,0.0);
    txphase(one);
    /* start of the CPMG period 1    */
    if(ncyc_cp == 1) {
        delay(tau_cpmg - (2.0/PI)*pw_cp);
        rgpulse(2.0*pw_cp,one,0.0,0.0);
        delay(tau_cpmg);
    }
    if(ncyc_cp == 2) {
        delay(tau_cpmg - (2.0/PI)*pw_cp);
        rgpulse(2.0*pw_cp,one,0.0,0.0);
        delay(tau_cpmg);
        delay(tau_cpmg);
        rgpulse(2.0*pw_cp,one,0.0,0.0);
        delay(tau_cpmg);
    }
    if(ncyc_cp > 2) {
        delay(tau_cpmg - (2.0/PI)*pw_cp);
        rgpulse(2.0*pw_cp,one,0.0,0.0);
        delay(tau_cpmg);
        initval(ncyc_cp-2,v4);
        loop(v4,v5);
          delay(tau_cpmg);
          rgpulse(2.0*pw_cp,one,0.0,0.0);
          delay(tau_cpmg);
         endloop(v5);
        delay(tau_cpmg);
        rgpulse(2.0*pw_cp,one,0.0,0.0);
        delay(tau_cpmg);
      }
    txphase(t4); decphase(zero);
    rgpulse(2.0*pw_cp,t4,2.0e-6,2.0e-6);
    txphase(one);
    if(ncyc_cp == 1) {
        delay(tau_cpmg);
        rgpulse(2.0*pw_cp,one,0.0,0.0); txphase(one);
        delay(tau_cpmg - 2.0/PI*pw_cp);
    }
    if(ncyc_cp == 2) {
        delay(tau_cpmg);
       rgpulse(2.0*pw_cp,one,0.0,0.0);
       delay(tau_cpmg);
       delay(tau_cpmg);
       rgpulse(2.0*pw_cp,one,0.0,0.0); txphase(one);
       delay(tau_cpmg - 2.0/PI*pw_cp);
   }
   if(ncyc_cp > 2) {
       delay(tau_cpmg);
       rgpulse(2.0*pw_cp,one,0.0,0.0);
       delay(tau_cpmg);
       initval(ncyc_cp-2,v4);
       loop(v4,v5);
         delay(tau_cpmg);
         rgpulse(2.0*pw_cp,one,0.0,0.0);
         delay(tau_cpmg);
        endloop(v5);
       delay(tau_cpmg);
       rgpulse(2.0*pw_cp,one,0.0,0.0); txphase(one);
       delay(tau_cpmg - 2.0/PI*pw_cp);
     }
     rgpulse(pw_cp,zero,0.0,0.0);
     delay(2.0e-6);
     rgradient('z',gzlvl5);
     delay(gt5);
     rgradient('z',0.0);
     delay(250.0e-6);
     obspower(tpwr);
     rgpulse(pw,zero,4.0e-6,0.0);
     decpower(d_reb);
   delay(2.0e-6);
   rgradient('z',gzlvl6);
   delay(gt6);
   rgradient('z',0.0);
   delay(150.0e-6);
   delay(taua - POWER_DELAY - gt6 - 152e-6
          - WFG2_START_DELAY - 0.5*pwc_reb);
   simshaped_pulse("hard",sh_reb,2.0*pw,pwc_reb,zero,zero,0.0,0.0);
   delay(taua - 0.5*pwc_reb - WFG2_STOP_DELAY - 2.0*POWER_DELAY - gt6 - 152e-6);
     rlpower(dpwr,DODEV); /* Set power for decoupling */
     rlpower(dpwr2,DO2DEV); /* Set power for decoupling */
     delay(2.0e-6);
     rgradient('z',gzlvl6);
     delay(gt6);
     rgradient('z',0.0);
     delay(150.0e-6);
     rgpulse(pw,zero,0.0,0.0);
/*     rcvron();  */          /* Turn on receiver to warm up before acq */
/* BEGIN ACQUISITION */
status(D);
     setreceiver(t5);
}
