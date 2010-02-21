#ifndef LINT
static char SCCSid[] = "@(#)DEPT.c 17.2 05/17/04 Copyright (c) 1991-1998 Varian Assoc.,Inc. All Rights Reserved";
#endif
/* 
 * Varian Assoc.,Inc. All Rights Reserved.
 * This software contains proprietary and confidential
 * information of Varian Assoc., Inc. and its contributors.
 * Use, disclosure and reproduction is prohibited without
 * prior consent.
 */
/* DEPT - distortionless enhancement by polarization transfer

	Parameters:
	   j1xh : 	CH coupling constant
	   mult :	Multiplication factor for the theta pulse
	   		 (set as an array = 0.5,1,1,1.5 for editing)
	   pplvl:	proton pulse level
	   pp	:	proton pulse width

KrishK	-	Last revision	: June 1997
*/


#include <standard.h>

pulsesequence()
{
   double          j1xh,
                   mult,
                   pp,
		   pplvl,
                   tau;


   j1xh      = getval("j1xh");
   mult   = getval("mult");
   pp     = getval("pp");
   pplvl = getval("pplvl");
   tau    = 1.0 / (2.0 * j1xh);

/* CHECK CONDITIONS */
   if ((rof1 + pp) < 2.0e-5)
      rof1 = (2.0e-5) - pp;
   if (rof1 < 1.0e-5)
      rof1 = 1.0e-5;
   if (rof2 == 0.0)
      rof2 = 1.0e-5;

   if ((dm[0] == 'y') || (dm[1] == 'y'))
   {
      fprintf(stdout, "Decoupler must be set as dm='nny'.\n");
      psg_abort(1);
   }
   if ((dmm[0] != 'c') || (dmm[1] != 'c'))
   {
      fprintf(stdout, "Decoupler must be set as dmm='ccf' or dmm='ccw'.\n");
      psg_abort(1);
   }


/* PHASECYCLE CALCULATIONS */
   hlv(ct, v5);
   hlv(v5, v4);
   hlv(v4, v4);
   hlv(v4, v7);
   hlv(v7, v2);
   mod2(v7, v7);
   dbl(ct, v1);
   assign(v1, v10);
   add(v1, v7, v1);
   dbl(v5, v9);
   add(v9, v10, v10);
   add(v5, v7, v5);
   add(one, v5, v6);
   mod2(v4, v4);
   dbl(v4, v4);
   add(v4, v7, v4);
   add(v4, v10, v10);
   dbl(v2, v9);
   add(v9, v10, v10);
   add(v2, v7, v2);
   add(one, v2, v3);
   mod4(v10, oph);
   add(one, v7, v8);


/* ACTUAL PULSESEQUENCE BEGINS */
   status(A);
      decpower(pplvl);
      delay(d1);

   status(B);
      decpulse(pp, v1);
      if (pw > 2.0*pp)
         delay(tau - 0.5*pp - 0.5*pw - rof2 - rof1);
      else
         delay(tau - 1.5*pp - rof2 - rof1);
      simpulse(pw, 2.0 * pp, v4, v3, rof1, rof2);

      if (mult * pp > 2.0 * pw)
         delay(tau - pp - pw - 0.5*mult*pp - rof2 - 2.1e-5);
      else
         delay(tau - pp - 2.0 * pw - rof2 - 2.1e-5);

      rgpulse(pw, v5, 2.0e-5, 0.0); 
      rcvroff();
      simpulse(2.0 * pw, mult * pp, v6, v8, 1.0e-6, 0.0);
      rgpulse(pw, v5, 1.0e-6, rof2);
      rcvron();

      decpower(dpwr);

      if (mult * pp > 2.0 * pw)
         delay(tau - pw - 0.5*mult*pp - 1.0e-6 - rof2);
      else
         delay(tau - 2.0*pw - 1.0e-6 - rof2);

   status(C);
      decphase(zero);
}
