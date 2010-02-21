/* this pulse sequence is provided on as-is basis
with no guarantees of fitness, merchantability, etc.
you are responsible for setting correct parameters
*/

#include <standard.h>

pulsesequence()
{

	double		gzlvl0,gt0,gzlvl1,gt1,gzlvl2,gt2,gzlvl3,gt3,gzlvl4,gt4; 
	double		gzlvl5,gt5,gzlvl6,gt6,gzlvl7,gt7,gzlvl8,gt8,gstab,gdelay;
	double		tCH,pwCBIP,pwC,pwHBIP,tpwrs,cpwr,npwr,pws;

	char    	shp_CBIP[MAXSTR],shp_HBIP[MAXSTR],shp[MAXSTR];

	int		iphase,icosel;	

	pwC    = getval("pwC");
	pwCBIP = getval("pwCBIP");
	pwHBIP = getval("pwHBIP");
        cpwr   = getval("cpwr");
        npwr   = getval("npwr");
	tpwrs  = getval("tpwrs");
	pws    = getval("pws");
	tCH    = getval("tCH");				/*tCH = 1/4*JCH */

        gzlvl0 = getval("gzlvl0");
	gt0    = getval("gt0");
	gzlvl1 = getval("gzlvl1");
        gt1    = getval("gt1");
	gzlvl2 = getval("gzlvl2");
        gt2    = getval("gt2");
	gzlvl3 = getval("gzlvl3");
        gt3    = getval("gt3");
	gzlvl4 = getval("gzlvl4");
        gt4    = getval("gt4");
	gzlvl5 = getval("gzlvl5");
        gt5    = getval("gt5");
        gzlvl6 = getval("gzlvl6");
        gt6    = getval("gt6");
        gzlvl7 = getval("gzlvl7");
        gt7    = getval("gt7");
        gzlvl8 = getval("gzlvl8");
        gt8    = getval("gt8");
	gstab  = getval("gstab");
	gdelay = getval("gdelay");

	getstr("shp_CBIP",shp_CBIP);
	getstr("shp",shp);
	getstr("shp_HBIP",shp_HBIP);	

	iphase  = (int)(getval("phase")+0.5);
	if (iphase == 1)
	      {
		assign(two,v4);
		icosel = -1;
	      }
	else	
	      {
		assign(zero,v4);
		icosel = 1;   
	      }
	
	initval(2.0*(double)(((int)(d2*getval("sw1")+0.5)%2)),v11);	

	mod2(ct,v3); 				/* v3  = 0101 0101 0101 ...*/
	dbl(v3,v3);           	         	/* oph = 0202 0202 0202 ...*/
	assign(v3,oph);

	mod4(ct,v1);				/* v1  = 0123 0123 0123 ... */ 

	assign(zero,v2);
	add(v11,v2,v2);
	add(v11,oph,oph);		

status(A);

	delay(d1/2.0);				/* equilibrium period */

	obspower(tpwr);       obsoffset(tof);
        decpower(cpwr);       decoffset(dof);
        dec2power(npwr);      dec2offset(dof2);

	delay(d1/2.0);

	rcvroff();

	/*decrgpulse(pwC,zero,rof1,rof2);
        zgradpulse(gzlvl0,gt0);
        delay(gstab);

        decrgpulse(pwC,zero,rof1,rof2);
        zgradpulse(gzlvl1,gt1);
        delay(gstab);*/

status(B);

	rgpulse(pw,zero,rof1,rof2);		/* Start of Pulse Sequence */
	shaped_pulse(shp_HBIP,pwHBIP,zero,rof1,rof2);
	zgradpulse(gzlvl2,gt2);
        delay(gstab);

	delay(tCH+(4.0*pw/3.1416)-rof1-rof2-gt2-gstab);

	simshaped_pulse(shp_HBIP,shp_CBIP,pwHBIP,pwCBIP,zero,v2,rof1,rof2);

        zgradpulse(gzlvl2,gt2);
        delay(gstab);    	

	delay(tCH-rof1-rof2-gt2-gstab);						

	rgpulse(pw,one,rof1,rof2);
	zgradpulse(gzlvl3,gt3);
	delay(gstab);

	decrgpulse(pwC,v2,rof1,rof2);		/* real Time */
	delay(d2/2.0);
	shaped_pulse(shp_HBIP,pwHBIP,zero,rof1,rof2);
	delay(d2/2.0);

	zgradpulse(icosel*gzlvl4,gt4);			/* CLUB*/
	delay(gdelay);
	decshaped_pulse(shp_CBIP,pwCBIP,zero,rof1,rof2);
	zgradpulse(-icosel*gzlvl4,gt4);
	delay(gdelay);
	delay(2.0*rof1+2.0*rof2+4.0*pwC/3.1416+pwHBIP+WFG_START_DELAY+WFG_STOP_DELAY);
	zgradpulse(-icosel*gzlvl5,gt5);
	delay(gdelay);
	decshaped_pulse(shp_CBIP,pwCBIP,v1,rof1,rof2);
	zgradpulse(icosel*gzlvl5,gt5);
	delay(gdelay);

	simpulse(pw,pwC,zero,v4,rof1,rof2);
	simshaped_pulse(shp_HBIP,shp_CBIP,pwHBIP,pwCBIP,zero,zero,rof1,rof2);
	zgradpulse(gzlvl6,gt6);
        delay(gstab);
	delay(tCH+(4.0*pw/3.1416)-gt6-gstab-rof1-rof2);
	simshaped_pulse(shp_HBIP,shp_CBIP,pwHBIP,pwCBIP,zero,zero,rof1,rof2);
	zgradpulse(gzlvl6,gt6);
        delay(gstab);
        delay(tCH-rof1-rof2-gt6-gstab);

	simpulse(pw,pwC,one,one,rof1,rof2);
	simshaped_pulse(shp_HBIP,shp_CBIP,pwHBIP,pwCBIP,zero,zero,rof1,rof2);
	zgradpulse(gzlvl7,gt7);
        delay(gstab);
        delay(tCH+(4.0*pw/3.1416)-gt7-gstab-rof1-rof2);
        simshaped_pulse(shp_HBIP,shp_CBIP,pwHBIP,pwCBIP,zero,zero,rof1,rof2);
        zgradpulse(gzlvl7,gt7);
        delay(gstab);
	delay(tCH-rof1-rof2-gt7-gstab);

	rgpulse(pw,zero,rof1,rof2);
        shaped_pulse(shp_HBIP,pwHBIP,zero,rof1,rof2);
	zgradpulse(-(gzlvl8/2.0),gt8);
	delay(gstab+(2.0*pw/3.1416));
	dec2power(dpwr2);
	shaped_pulse(shp_HBIP,pwHBIP,zero,rof1,rof2);
	zgradpulse(gzlvl8/2.0,gt8);
	delay(gstab);
	decpower(dpwr);

	rcvron();
	
status(C);	
}
