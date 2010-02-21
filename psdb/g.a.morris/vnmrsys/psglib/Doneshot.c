#ifndef LINT
static char SCCSid[] = "@(#)Doneshot.c 18.2 12/23/04 Copyright (c) 1991-2004 Varian Inc. All Rights Reserved";
#endif
/* 
 * Varian Inc. All Rights Reserved.
 * This software contains proprietary and confidential
 * information of Varian Inc. and its contributors.
 * Use, disclosure and reproduction is prohibited without
 * prior consent.
 */
/**********************************************************
	Doneshot - Oneshot DOSY sequence based on modifications to the
	Bipolar pulse stimulated echo (BPPSTE) sequence

Parameters:
        delflag   - 'y' runs the Doneshot sequence
                    'n' runs the normal s2pul sequence
	del       -  the actual diffusion delay
	gt1       - total diffusion-encoding pulse width
        gzlvl1    - diffusion-encoding pulse strength
        gstab     - gradient stabilization delay (~0.0002-0.0003 sec)
        gt3       - spoiling gradient duration (in sec)
        gzlvl3    - spoiling gradient strength (destroys transverse 
			magnetisation during the diffusion delay)
        gzlvl_max - maximum accepted gradient strength
                       32767 with PerformaII, 2047 with PerformaI
        kappa     - unbalancing factor between bipolar pulses as a
                       proportion of gradient strength (~0.2)
	fn2D      - Fourier number to up the 2D display in F2
        wet       - 'y' turns wet flag on
        satmode   - 'y' turns on presaturation


The parameters for the heating gradients (gt4, gzlvl4) are calculated
in the sequence. They cannot be set directly.
tau defined as time between the mid-points of the bipolar diffusion encoding
gradient pulses

	Constant energy dissipation calculation corrected 8vii99

17xi08 GAM	Add warning if phasecycleflag='n'  !
8iv09 GAM	Just remove phasecycleflag, no longer needed
8iv09 GAM	Correct value for Dtau (remove spurious gt1/2)


**********************************************************/

#include <standard.h>

pulsesequence()
{
double	kappa     = getval("kappa"), 
	gzlvl1    = getval("gzlvl1"),
	gzlvl3    = getval("gzlvl3"),
        gzlvl_max = getval("gzlvl_max"),
	gt1       = getval("gt1"),
	gt3       = getval("gt3"),
	del       = getval("del"),
        gstab     = getval("gstab"),
	gzlvl4,gt4,Dtau,Ddelta,dosytimecubed, dosyfrq,
        satpwr = getval("satpwr"),
        satdly = getval("satdly"),
        satfrq = getval("satfrq");
char delflag[MAXSTR], wet[MAXSTR],satmode[MAXSTR];

   gt4 = 2.0*gt1;
   getstr("delflag",delflag);
   getstr("wet",wet);
   getstr("satmode",satmode);

/* Decrement gzlvl4 as gzlvl1 is incremented, to ensure constant 
   energy dissipation in the gradient coil 
   Current through the gradient coil is proportional to gzlvl */

   gzlvl4 = sqrt(2.0*gt1*(1+3.0*kappa*kappa)/gt4) * 
            sqrt(gzlvl_max*gzlvl_max/((1+kappa)*(1+kappa))-gzlvl1*gzlvl1);

/* In pulse sequence, del>4.0*pw+3*rof1+2.0*gt1+5.0*gstab+gt3 */

   if ((del-(4*pw+3.0*rof1+2.0*gt1+5.0*gstab+gt3)) < 0.0)
   {  del=(4*pw+3.0*rof1+2.0*gt1+5.0*gstab+gt3);
      text_message("Warning: del too short; reset to minimum value");
   }

   if ((d1 - (gt3+gstab) -2.0*(gt4/2.0+gstab)) < 0.0)
   {  d1 = (gt3+gstab) -2.0*(gt4/2.0+gstab);
      text_message("Warning: d1 too short;  reset to minimum value");
   }

   if ((gzlvl1*(1+kappa)) > gzlvl_max)
   {  abort_message("Max. grad. amplitude exceeded: reduce either gzlvl1 or kappa");
   }

   if (ni > 0.0)
   {  abort_message("This is a 2D, not a 3D dosy sequence: please set ni to zero");
   }

   Ddelta=gt1;
   Dtau=2.0*pw+gstab+rof1;
   dosyfrq = sfrq;
   dosytimecubed = Ddelta*Ddelta*(del+(Ddelta/6.0) *
                   (kappa*kappa-2.0) +
                   (Dtau/2.0)*(kappa*kappa-1.0));
   putCmd("makedosyparams(%e,%e)\n",dosytimecubed,dosyfrq);

/* phase cycling calculation */

   if(delflag[0]=='y')
   {
	 mod2(ct,v1); dbl(v1,v1);  hlv(ct,v2);
	 mod2(v2,v3); dbl(v3,v3);  hlv(v2,v2);
	 mod2(v2,v4); add(v1,v4,v1);			/*    v1      */
	 hlv(v2,v2);  add(v2,v3,v4);			/*    v4      */
	 hlv(v2,v2);  mod2(v2,v3); dbl(v3,v5);
	 hlv(v2,v2);  mod2(v2,v3); dbl(v3,v3);		/*    v3      */
	 hlv(v2,v2);  mod2(v2,v6); add(v5,v6,v5);	/*    v5      */
	 hlv(v2,v2);  mod2(v2,v2); dbl(v2,v2);		/*    v2      */
	 assign(v1,oph);  dbl(v2,v6);      sub(oph,v6,oph);
	 add(v3,oph,oph); sub(oph,v4,oph); dbl(v5,v6);
	 add(v6,oph,oph); mod4(oph,oph);                /* receiver phase */
   }
   status(A);

   if(delflag[0]=='y')
   {  zgradpulse(-1.0*gzlvl4,gt4/2.0);	/* 1st dummy heating pulse */
      delay(gstab);

      zgradpulse(gzlvl4,gt4/2.0);	/* 2nd dummy heating pulse */
      delay(gstab);

      delay(d1 - (gt3+gstab) -2.0*(gt4/2.0+gstab));

      zgradpulse(-1.0*gzlvl3,gt3);	/* Spoiler gradient balancing pulse */
      delay(gstab); 
   }
   else
   if (satmode[0] == 'y')
       {
       if (d1 - satdly > 0)
         delay(d1 - satdly);
       else
       delay(0.02);
       obspower(satpwr);
       txphase(v1);
        if (satfrq != tof)
         obsoffset(satfrq);
        rgpulse(satdly,zero,rof1,rof1);
        if (satfrq != tof)
         obsoffset(tof);
       obspower(tpwr);
       delay(1.0e-5);
      }
     else
     {  delay(d1); }

   if (wet[0] == 'y')     wet4(zero,one);


   status(B); /* first part of sequence */
   if (delflag[0]=='y')
   {
      if (gt1>0 && gzlvl1>0)
      {
   	 rgpulse(pw, v1, rof1, 0.0);		/* first 90, v1 */

	 zgradpulse(gzlvl1*(1.0-kappa),gt1/2.0); /*1st main gradient pulse*/
   	 delay(gstab);
	 rgpulse(pw*2.0, v2, rof1, 0.0);	/* first 180, v2 */

	 zgradpulse(-1.0*gzlvl1*(1.0+kappa),gt1/2.0); /*2nd main grad. pulse*/
   	 delay(gstab);
   	 rgpulse(pw, v3, rof1, 0.0);		/* second 90, v3 */

	 zgradpulse(gzlvl1*2.0*kappa,gt1/2.0);  /* Lock refocussing pulse*/
   	 delay(gstab);

	 zgradpulse(gzlvl3,gt3); /* Spoiler gradient balancing pulse */
   	 delay(gstab);

	 delay(del-4.0*pw-3.0*rof1-2.0*gt1-5.0*gstab-gt3); /* diffusion delay */

	 zgradpulse(2.0*kappa*gzlvl1,gt1/2.0);	/*Lock refocussing pulse*/
   	 delay(gstab);
   	 rgpulse(pw, v4, rof1, 0.0);		/* third 90, v4 */

	 zgradpulse(gzlvl1*(1.0-kappa),gt1/2.0); /*3rd main gradient pulse*/
   	 delay(gstab);
	 rgpulse(pw*2.0, v5, rof1, rof2);	/* second 180, v5 */

	 zgradpulse(-1.0*(1.0+kappa)*gzlvl1,gt1/2.0); /*4th main grad. pulse*/
   	 delay(gstab);
      }
   }
   else
      rgpulse(pw,oph,rof1,rof2);
   status(C);
}

/****************************************************************************
256 (or 16,32,64,128) steps phase cycle
Order of cycling: v1 (0,2), v4(0,2), v1(1,3), v4(1,3), v5(0,2), v3(0,2),
v5(1,3) v2(0,2)

receiver = v1-2*v2+v3-v4+2*v5

Coherence pathway :

	90	180	 90		90	180	  Acq.
+2
                  ________                _______
+1               |        |              |       |
                 |        |              |       |
 0------|        |        |--------------|       |
        |        |                               |
-1      |________|                               |____________________

-2
	v1       v2        v3             v4      v5
****************************************************************************/
