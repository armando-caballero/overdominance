/* purgingf13class.c (1-sod1 1 1-sod2)*/

/* ***************************************************** */

#include "libhdr"
#include "ranlib.h"
#define NN 1001  /* max number of NIND and gen */
#define MM 2001  /* max number of NCRO */
#define maxmpt 5001
#define normalthreshold 30

int muts, mutsL, mutsOD, lastinmutantspoissontable, lastinmutantspoissontableL, lastinmutantspoissontableOD, lastinrecpoissontable;
int marker, NIND, NCRO, NLOCI, TOTLOCI, numSegSD;
int gen, generations, i, k, l, m, lin, rep, lines, replicates, dom;
int mutants_ocurred, RM[NN], initialgen[MM][31], p1, p2;
int ran_i, ran_k, ran_l, ran_h, countNoSS, NoSS_k[15001], NoSS_l[15001];
int gm[NN][MM][2], sm[NN][MM][2];
int NINDNP, gmNP[5001][MM][2];
int initialgenNP[MM][31];
int rnd;

double addedfixed_sNP, sNP[MM][31], hsNP[MM][31];
double genvalm_s[NN], pm_s[NN], numAaSD;
double genvalm_s_F[NN], pm_s_F[NN];
double va, vd, id, d_a, alfa_a;
double va_let, vd_let, id_let, d_a_let, alfa_a_let;
double va_od, vd_od, id_od, d_a_od, alfa_a_od;
double va_del, vd_del, id_del, d_a_del, alfa_a_del;
double L, PS, Lambda_s, Lambda_L, Lambda_OD;
double AA, Aa, aa, q[MM][31];
double ave_s, addedfixed_s, alpha_s, beta_s, k_s, ave_hs, zs;
double mutantspoissontable[maxmpt], mutantspoissontableL[maxmpt], mutantspoissontableOD[maxmpt];
double s[MM][31], hs[MM][31], qbar[NN][MM];
double H, Q, QS[NN], HO[NN], HS[NN], Fstjunk, HTjunk;
double recpoissontable[maxmpt];
double sod1, sod2;
double fyang2[NN];

double s_real[MM][31], hs_real[MM][31];
double s1_real[MM][31], s2_real[MM][31];

struct acc SK2[NN], gmean_s[NN], gmean_s_F[NN], Wmax[NN], gvar_s[NN], AVE_VA[NN], AVE_VD[NN], AVE_ID[NN];
struct acc AVE_VA_let[NN], AVE_VD_let[NN], AVE_ID_let[NN], AVE_VA_od[NN], AVE_VD_od[NN], AVE_ID_od[NN], AVE_VA_del[NN], AVE_VD_del[NN], AVE_ID_del[NN];
struct acc del[NN], let[NN], od[NN], qdel[NN], qlet[NN], qod[NN];
struct acc QSrep[NN], HOrep[NN], HSrep[NN], HTrep[NN];
struct acc FST[NN], NeH[NN], NeFst[NN];
struct acc numberofrecom;
struct acc ID_HOM[NN];
struct acc fix_sd[NN], lost_sd[NN]; 

struct acc sdel_X[NN], hsdel_X[NN], sdel_real[NN], hsdel_real[NN], sod1_real[NN], sod2_real[NN];
struct acc AVE_VA_del_real[NN], AVE_VD_del_real[NN], AVE_ID_del_X[NN], AVE_ID_del_real[NN], AVE_VA_od_real[NN], AVE_VD_od_real[NN], AVE_ID_od_X[NN], AVE_ID_od_real[NN];

// ADDITIONAL EXPERIMENT
int I, J;
int gmG0[NN][MM][2], gmG1[NN][MM][2], gmG2[NN][MM][2]; 
double pm_sG0[NN], pm_sG1[NN], pm_sG2[NN], mean_parents_G0[NN], mean_family_G1[NN], mean_family_G2[NN]; 
struct acc AVE_VAeff[NN], AVE_VDeff[NN], AVE_B_FS[NN];
struct acc IDFyang2[NN];

void natural_population(), sample(), phenotypeB();
void recombination_masks(), genotypic_values(), real_values();
void neutral_genes(), selected_genes();
void mutation_neutral(), mutation_selected(), disorder_NoSS();
void mating();

FILE *fptr, *fgen, *frep, *fdat, *fpop, *fopen();


/* ***************************************************** */


main()
{
	fptr = fopen ("dfilename.dat","w");
	fgen = fopen ("genfile.dat","w");

	getinputs();
	headings();
	recombination_masks(RM);
	initialize();
	natural_population(gmNP, sNP, hsNP, initialgenNP, &NINDNP);

	for (rep=1; rep<=replicates; rep++)
	{
		frep = fopen ("repfile.dat","a");
		fprintf (frep,"replicate %d\n", rep);
		fclose(frep);

		//if (evel!=0) fprintf (fptr,"replicate %d\n", rep);

		for (gen=0; gen<=generations; gen++)
		{
			QS[gen]=0.0;
			HO[gen]=0.0;
			HS[gen]=0.0;
			for (k=0; k<NCRO; k++)	qbar[gen][k]=0.0;
		}

		for (lin=1; lin<=lines; lin++)
		{
			natural_data();
			sample(gmNP, gm);

			for (gen=0; gen<=generations; gen++)
			{
				neutral_genes (gm, q, initialgen);
				selected_genes (gm, q, s, hs, initialgen);
				mutation_neutral (gm, q, initialgen);
				mutation_selected (gm, q, s, hs, initialgen);
				//if (tracelevel!=0) dumpoffspringaftermutation();
				genotypic_values (gm, s, hs, genvalm_s, RM);
				phenotypeB (genvalm_s, pm_s);
				real_values(s, hs, genvalm_s);
 
				// ADDITIONAL EXPERIMENT
				generation1();
				generation2();
                
				//if (tracelevel!=0) dumpphenotypes();
				mating (gm, sm, pm_s);
				//if (tracelevel!=0) dumpoffspring();
			}
		}
		for (gen=0; gen<=generations; gen++)
		{
			accum(&HOrep[gen], HO[gen]/lines);
			accum(&QSrep[gen], QS[gen]/lines);
			HS[gen] = HS[gen]/lines ;
			accum(&HSrep[gen], HS[gen]);
			accum(&NeH[gen], 1.0 / ( 2.0 * (1.0 - pow(HS[gen] / HS[0], 1.0/gen))));

			HTjunk = 0.0;
			for (k=0; k<NCRO; k++)	HTjunk += 2.0 * (qbar[gen][k]/lines) * ( 1.0 - (qbar[gen][k]/lines));
			// Nei and Chesser correction
			HTjunk = HTjunk/NCRO + (HS[gen]/(NIND*lines)) - (HS[gen]/(2.0*NIND*lines)) ;
			accum(&HTrep[gen], HTjunk);
			Fstjunk = (HTjunk - HS[gen]) / HTjunk;
			accum(&FST[gen], Fstjunk);
			accum(&NeFst[gen], 1.0 / ( 2.0 * (1.0 - pow(1.0-Fstjunk, 1.0/(double)gen))));
		}
  	}
	printout();

	writeseed();
}


/* ***************************************************** */

getinputs()
{
	tracestart();
	getseed();
	getintandskip("NIND (max 1000):",&NIND,2,1000);
	getrealandskip("probability of selfing (Random:99):",&PS,0.0,99.0);
	getintandskip("J family size (max 10):",&J,1,10);
	getrealandskip("Length of genome in Morgans (99:FreeRecom) :",&L,0.0,99.0);
	getintandskip("NCRO (min 1, max 2000):",&NCRO,1,2000);
	getintandskip("NLOCI (first is neutral)(min 2, max 30):",&NLOCI,2,30);
	TOTLOCI = NCRO * NLOCI;
	getrealandskip("Lambda_s :",&Lambda_s,0.0,(double)infinity);
	getrealandskip("Lambda_L (s=1,h=0.02):",&Lambda_L,0.0,(double)infinity);
	getrealandskip("Lambda_OD :",&Lambda_OD,0.0,(double)infinity);
	getrealandskip("sod1 :",&sod1,0.0,1.0);
	getrealandskip("sod2 :",&sod2,0.0,1.0);

	getrealandskip("Beta_s :",&beta_s,0.0,(double)infinity);

	getrealandskip("Average |s| :",&ave_s,0.0,1.0);
	alpha_s = beta_s / ave_s;

//	getintandskip("dom. model (0:constant 1:Deng 2:CK94(gamma):", &dom,0,2);
//	getrealandskip("hs (mod 0) or k_s (mod 1):", &k_s,0.0,(double)infinity);
	dom = 2;
	getrealandskip("Average h (mod 2):",&ave_hs,0.0,(double)infinity);
	if (dom==2)	k_s = alpha_s * (pow((2.0*ave_hs),((-1.0)/beta_s))-1.0);

	getintandskip("Number of generations :",&generations,2,1000);
	getintandskip("Number of lines :",&lines,1,infinity);
	getintandskip("Number of replicates :",&replicates,1,infinity);
}


/* **************************************************** */


headings()
{
	fgen = fopen ("genfile.dat","a");
	fprintf(fgen,"\nN=%d  PS=%f N.S-LOCI=%d    N.N-LOCI=%d    lines=%d   reps=%d   gens=%d\n    L=%f   Lambda_s=%6.4f   Lambda_L=%6.4f   Lambda_OD=%6.5f   sod1=%6.5f   sod2=%6.5f\n    ave_s=%6.4f   beta_s=%6.4f   alpha_s=%6.4f  k_s=%6.4f   ave_hs=%6.4f\n", NIND, PS, TOTLOCI-NCRO, NCRO, lines, replicates, generations, L, Lambda_s, Lambda_L, Lambda_OD, sod1, sod2, ave_s, beta_s, alpha_s, k_s, ave_hs);
	fprintf(fgen,"gen	W	    W_hom	ID_hom	ID_yang Wmax	L	    V(W)	Sk2	    VA	    VAeff	VAdel	VAlet	VAod    VD	    VDeff   VDdel    VDlet  VDod    B	    B_FS   Bdel	Blet	Bod	    FST	    NeVk	NeH	   NeFst\n");
	fclose(fgen);
}


/* ***************************************************** */


void recombination_masks (RM)
int RM[];
{
	for (l=0; l<NLOCI; l++)   RM[l]=pow(2.0,(double)l);
}


/* ***************************************************** */


initialize ()
{
    /* SELECTED LOCI WITH POISSON (2NL) NEW MUTATIONS */

    if ( (exp(-2.0*(double)NIND*Lambda_s) != 0.0)&&(2.0*(double)NIND*Lambda_s < normalthreshold) )
    generatepoissontable(2.0*(double)NIND*Lambda_s, &lastinmutantspoissontable, mutantspoissontable, maxmpt-1);

    if ( (exp(-2.0*(double)NIND*Lambda_L) != 0.0)&&(2.0*(double)NIND*Lambda_L < normalthreshold) )
    generatepoissontable(2.0*(double)NIND*Lambda_L, &lastinmutantspoissontableL, mutantspoissontableL, maxmpt-1);
    
    if ( (exp(-2.0*(double)NIND*Lambda_OD) != 0.0)&&(2.0*(double)NIND*Lambda_OD < normalthreshold) )
    generatepoissontable(2.0*(double)NIND*Lambda_OD, &lastinmutantspoissontableOD, mutantspoissontableOD, maxmpt-1);

	/* POISSON RECOMBINATION NUMBER */
	if ( (exp(-(double)L) != 0.0)&&((double)L < normalthreshold) )
	generatepoissontable((double)L, &lastinrecpoissontable, recpoissontable, maxmpt-1);
}


/* ***************************************************** */


void natural_population (int gmNP[][MM][2], double sNP[][31], double hsNP[][31], int initialgenNP[][31], int *NINDNP)
{
	int g0, g1, dinitialgen, dN;

	double dadds, dadda, ds, da, dhs, dha;

	/* ***** take genotypic values of natural population ***** */

	fpop=fopen("popfile","r");

	fscanf(fpop,"%d", &dN);
	*NINDNP = dN;

	for (i=0; i<dN; i++)
	for (k=0; k<NCRO; k++)
	{
		fscanf(fpop,"%d%d", &g0, &g1);
		gmNP[i][k][0] = g0;
		gmNP[i][k][1] = g1;
	}

	fclose(fpop);

	/* ***** take effects of genes ***** */

	fdat=fopen("datafile","r");

	fscanf(fdat,"%lf%lf", &dadds, &dadda);
	addedfixed_sNP=dadds;

	for (k=0; k<NCRO; k++)
	for (l=0; l<NLOCI; l++)
	{
		fscanf(fdat,"%lf%lf%lf%lf%d", &ds, &da, &dhs, &dha, &dinitialgen);
		sNP[k][l] = ds;
		hsNP[k][l] = dhs;

		initialgenNP[k][l] = dinitialgen;

//		if (tracelevel!=0)    if ((k<5)&&(l<5)) fprintf(fptr,"\n k=%d l=%d sNP=%f hsNP=%f initialgenNP=%d", k, l, sNP[k][l], hsNP[k][l], initialgenNP[k][l]);

		if (sNP[k][l]>0.99999) sNP[k][l]=0.99999;
		if (sNP[k][l]<(-1.0)) sNP[k][l]=(-1.0);
	}

	fclose(fdat);
}


/* ***************************************************** */


natural_data ()
{
	addedfixed_s = addedfixed_sNP;

	for (k=0; k<NCRO; k++)
	for (l=0; l<NLOCI; l++)
	{
		s[k][l] = sNP[k][l];
		hs[k][l] = hsNP[k][l];

		initialgen[k][l] = initialgenNP[k][l];
	}
}


/* ***************************************************** */
void sample (int gmNP[][MM][2], int gm[][MM][2])
{
	int g;

	for (i=0; i<NINDNP; i++)
	{
		ran_i = (int)(uniform()*NINDNP);
		for (k=0; k<NCRO; k++)
		{
			g=gmNP[i][k][0]; gmNP[i][k][0]=gmNP[ran_i][k][0]; gmNP[ran_i][k][0]=g;
			g=gmNP[i][k][1]; gmNP[i][k][1]=gmNP[ran_i][k][1]; gmNP[ran_i][k][1]=g;
		}
	}

	for (i=0; i<NIND; i++)
	for (k=0; k<NCRO; k++)
	{
		gm[i][k][0]=gmNP[i][k][0];
		gm[i][k][1]=gmNP[i][k][1];
	}
}

/* ***************************************************** */


void neutral_genes (int gm[][MM][2], double q[][31], int initialgen[][31])
{
	double Ho;

	Ho=0.0; H=0.0; Q=0.0;

	for (k=0; k<NCRO; k++)
	if(initialgen[k][0] != (-99))
	{
	    AA=0.0; Aa=0.0; aa=0.0;

	    for (i=0; i<NIND; i++)
	    {
	        if (((gm[i][k][0] & RM[0])==RM[0])&&((gm[i][k][1] & RM[0])==RM[0]))	aa+=1.0;
	        else    if (((gm[i][k][0] & RM[0])!=RM[0])&&((gm[i][k][1] & RM[0])!=RM[0])) AA+=1.0;
	        else		Aa+=1.0;
	    }

		Ho += Aa/(double)NIND;
		q[k][0] = (aa/(double)NIND)+(Aa/(2.0*(double)NIND));
		H += 2.0 * q[k][0] * (1.0 - q[k][0]);
		Q += q[k][0];
		qbar[gen][k] += q[k][0];

//	     if (tracelevel!=0)    fprintf(fptr,"%d\t0\t%1.0f\t%1.0f\t%1.0f\t%f\t%d\n",k,AA,Aa,aa,q[k][0],initialgen[k][0]);
	}
	HO[gen] += Ho/NCRO;
	QS[gen] += Q/NCRO;
	// Nei and Chesser correction
	HS[gen] += (NIND/(NIND-1.0)) * ((H/NCRO) - ((H/NCRO)/(2.0*NIND))) ;
}


/* ***************************************************** */


void selected_genes (int gm[][MM][2], double q[][31], double s[][31], double hs[][31], int initialgen[][31])
{
	for (k=0; k<NCRO; k++)
	for (l=1; l<NLOCI; l++)
	if (initialgen[k][l] != (-99))
	{
		AA=0.0; Aa=0.0; aa=0.0;

		for (i=0; i<NIND; i++)
		{
			if (((gm[i][k][0] & RM[l])==RM[l])&&((gm[i][k][1] & RM[l])==RM[l]))	aa+=1.0;
			else if (((gm[i][k][0] & RM[l])!=RM[l])&&((gm[i][k][1] & RM[l])!=RM[l]))	AA+=1.0;
	     	else	Aa+=1.0;
		}

		q[k][l] = (aa/(double)NIND)+(Aa/(2.0*(double)NIND));

		//if ((tracelevel!=0)&&(s[k][l]==(-1.0)))    fprintf(fptr,"\n q[%d][%d]=%f", k, l, q[k][l]);

		if ((q[k][l] > 0.0) && (q[k][l] < 1.0))
		{
			if (s[k][l]<= (-0.9))
			{
				accum (&let[gen], 1.0);
				accum (&qlet[gen], q[k][l]);
			}
			else if (hs[k][l]==99)
			{
				accum (&od[gen], 1.0);
				accum (&qod[gen], q[k][l]);
			}
			else
			{
				accum (&del[gen], 1.0);
				accum (&qdel[gen], q[k][l]);
			}
		}

		if (initialgen[k][l]==(-99)) /* gene non segregating */;
		else if (q[k][l]==0.0)
		{
			if (hs[k][l] == 99)
			{
				addedfixed_s *= (1.0 - sod1);
				accum (&lost_sd[gen], 1.0);
			}

			initialgen[k][l]=(-99);

			q[k][l]=0.0;
			s[k][l]=0.0; hs[k][l]=0.0;
		}
		else if (q[k][l]==1.0)
		{
			/* make it wild type */
			for (i=0; i<NIND; i++)
			{
				gm[i][k][0] = ( gm[i][k][0] & (~RM[l]) );
				gm[i][k][1] = ( gm[i][k][1] & (~RM[l]) );
			}
			initialgen[k][l]=(-99);
			if (hs[k][l] != 99)		addedfixed_s *= (1.0 + s[k][l]);
			if (hs[k][l] == 99)
			{
				addedfixed_s *= (1.0 - sod2);
				accum (&fix_sd[gen], 1.0);
			}

			q[k][l]=0.0;
			s[k][l]=0.0; hs[k][l]=0.0;
		}
	}
}


/* ***************************************************** */


void mutation_neutral (gm,q,initialgen)
int gm[][MM][2], initialgen[][31];
double q[][31];
{
    /* NEUTRAL GENES: (POISSON) 2N(Lambda_a) NEW MUTATIONS (TWO DIRECTIONS) */

    //if (tracelevel!=0)    fprintf(fptr,"\n New neutral mutants\n");

    muts = mutationnumber();

    for (m=0; m<muts; m++)
    {
		ran_i = (int)(uniform()*NIND);
		ran_k = (int)(uniform()*NCRO);
		ran_h = (int)(uniform()*2.0);
		if ( (gm[ran_i][ran_k][ran_h] & RM[0])==RM[0] )		gm[ran_i][ran_k][ran_h]=(gm[ran_i][ran_k][ran_h] & (~RM[0]));
		else		gm[ran_i][ran_k][ran_h]=(gm[ran_i][ran_k][ran_h] | RM[0]);
    }
}


/* ***************************************************** */


void mutation_selected (gm,q,s,hs,initialgen)
int gm[][MM][2], initialgen[][31];
double q[][31], s[][31], hs[][31];
{
    /* SELECTED GENES: (POISSON) 2N(Lambda_s + Lambda_L) NEW MUTATIONS */
   
    muts = mutationnumber();
    mutsL = mutationnumberL();
    mutsOD = mutationnumberOD();

    if ( (muts + mutsL + mutsOD) == 0 ) goto label;

	//if (tracelevel!=0)    fprintf(fptr,"\n New selected mutants  muts=%d  mutsL=%d\n", muts, mutsL);

    countNoSS = 0;

    for (k=0; k<NCRO; k++)
    {
		for (l=1; l<NLOCI; l++)
		{
			if (q[k][l]==0.0)
	   		{
				countNoSS += 1;
				NoSS_k[countNoSS-1] = k;
				NoSS_l[countNoSS-1] = l;
   	   		}
		}
    }

    //if (tracelevel!=0)    fprintf(fptr,"\n countNoSS=%d\n", countNoSS);

    mutants_ocurred = 0;

    if (countNoSS != 0)
    {
	disorder_NoSS (NoSS_k,NoSS_l);

	for (m=0; m<countNoSS; m++)
	{
		if (mutants_ocurred==(muts+mutsL+mutsOD))    goto label;

		ran_i = (int)(uniform()*NIND);
		ran_h = (int)(uniform()*2.0);
		gm[ran_i][NoSS_k[m]][ran_h]=(gm[ran_i][NoSS_k[m]][ran_h] | RM[NoSS_l[m]]);
	  	mutants_ocurred += 1;
		initialgen[NoSS_k[m]][NoSS_l[m]] = gen;

			/* ****** Lethal and OD mutations ****** */

			if(mutants_ocurred <= mutsL)
			{
				s[NoSS_k[m]][NoSS_l[m]] = (-1.0);
				hs[NoSS_k[m]][NoSS_l[m]] = 0.02;
			}
			else if(mutants_ocurred <= (mutsL+mutsOD))
			{
				hs[NoSS_k[m]][NoSS_l[m]] = 99;
			}
			else
			{
		     	/* ****** values of s, hs ****** */
				zs = gengam (alpha_s, beta_s);
				if(zs > 1.0)	zs = 1.0;
				s[NoSS_k[m]][NoSS_l[m]] = (-zs);
				if (s[NoSS_k[m]][NoSS_l[m]] >= (-0.42))	hs[NoSS_k[m]][NoSS_l[m]] = uniform() * exp(k_s*s[NoSS_k[m]][NoSS_l[m]]);
				else								hs[NoSS_k[m]][NoSS_l[m]] = uniform() * 0.04;
			}
			//if (tracelevel!=0)    fprintf(fptr,"New selected mutants  muts=%d  mutsL=%d  mutsOD=%d\n", muts, mutsL, mutsOD);
			//if (tracelevel!=0)    fprintf(fptr,"s=%f\n", s[NoSS_k[m]][NoSS_l[m]]);
			//if (tracelevel!=0)    fprintf(fptr,"ran_i=%d  k=%d  l=%d   ran_h=%d\n", ran_i, NoSS_k[m],  NoSS_l[m], ran_h);
		}
    }

    for (m=mutants_ocurred; m<muts; m++)
    {
		ran_i = (int)(uniform()*NIND);
		ran_k = (int)(uniform()*NCRO);
		do {ran_l = (int)(uniform()*NLOCI);}   while (ran_l== 0);
		ran_h = (int)(uniform()*2.0);

		gm[ran_i][ran_k][ran_h]=(gm[ran_i][ran_k][ran_h] | RM[ran_l]);

		//if (tracelevel!=0)    fprintf(fptr,"ran_i=%d  ran_k=%d  ran_l=%d  ran_h=%d\n", ran_i, ran_k, ran_l, ran_h);
    }
    label: /* end of mutations */;
}


/* ***************************************************** */


int mutationnumber ()
{
	int r;
	if ((2.0*(double)NIND*Lambda_s < normalthreshold) && (exp(-2.0*(double)NIND*Lambda_s) != 0.0) )
	{
		r= poisson(lastinmutantspoissontable, mutantspoissontable);
	}
	else r = (int)( normal(2.0*(double)NIND*Lambda_s, sqrt(2.0*(double)NIND*Lambda_s)) );
	return(r);
}


/* ***************************************************** */


int mutationnumberL ()
{
	int r;
	if ((2.0*(double)NIND*Lambda_L < normalthreshold) && (exp(-2.0*(double)NIND*Lambda_L) != 0.0) )
	{
		r= poisson(lastinmutantspoissontableL, mutantspoissontableL);
	}
	else r = (int)( normal(2.0*(double)NIND*Lambda_L, sqrt(2.0*(double)NIND*Lambda_L)) );
	return(r);
}


/* ***************************************************** */

int mutationnumberOD ()
{
	int r;
	if ((2.0*(double)NIND*Lambda_OD < normalthreshold) && (exp(-2.0*(double)NIND*Lambda_OD) != 0.0) )
	{
		r= poisson(lastinmutantspoissontableOD, mutantspoissontableOD);
	}
	else r = (int)( normal(2.0*(double)NIND*Lambda_OD, sqrt(2.0*(double)NIND*Lambda_OD)) );
	return(r);
}


/* ***************************************************** */


void disorder_NoSS (NoSS_k,NoSS_l)
int NoSS_k[], NoSS_l[];
{
	int a, b, rnd;

	for (i=0; i<countNoSS-1; i++)
	{
	   rnd=(int)(uniform()*(countNoSS-i));
	   a=NoSS_k[countNoSS-1-i]; NoSS_k[countNoSS-1-i]=NoSS_k[rnd]; NoSS_k[rnd]=a;
	   b=NoSS_l[countNoSS-1-i]; NoSS_l[countNoSS-1-i]=NoSS_l[rnd]; NoSS_l[rnd]=b;
	}
}


/* ***************************************************** */


dumpoffspringaftermutation()
{
	//if (tracelevel==0)   return (0);

	//fprintf(fptr,"\n Offspring after mutation (gm0 gm1)\n");
	//for (i=0; i<NIND; i++)   fprintf(fptr,"%d  %d\n",gm[i][0][0],gm[i][0][1]);
}


/* ***************************************************** */


void genotypic_values (int gm[][MM][2], double s[][31], double hs[][31], double genvalm_s[], int RM[])
{
	//if (tracelevel!=0)	fprintf(fptr,"\n Genotypic values \n");

    for (i=0; i<NIND; i++)
	{
		numAaSD = 0.0;

		genvalm_s[i] = addedfixed_s;
		genvalm_s_F[i] = addedfixed_s;
        
		if (uniform() < 0.5) rnd = 1;
		else				 rnd = 0;
        
		for (k=0; k<NCRO; k++)
		for (l=1; l<NLOCI; l++)
		if ((initialgen[k][l] != (-99)) && (hs[k][l] != 99))
		{
			if (((gm[i][k][0] & RM[l])==RM[l])&&((gm[i][k][1] & RM[l])==RM[l]))
			{
				genvalm_s[i] *= (1.0 + s[k][l]);
				genvalm_s_F[i] *= (1.0 + s[k][l]);
			}
			else if (((gm[i][k][0] & RM[l])!=RM[l])&&((gm[i][k][1] & RM[l])!=RM[l])) /* AA */;
			else if (((gm[i][k][0] & RM[l])==RM[l])&&((gm[i][k][1] & RM[l])!=RM[l]))
			{
				genvalm_s[i] *= (1.0 + (s[k][l]*hs[k][l]));
				if (rnd == 1)
				{
					genvalm_s_F[i] *= (1.0 + s[k][l]);
				}
				else		/*11*/;
			}
			else if (((gm[i][k][0] & RM[l])!=RM[l])&&((gm[i][k][1] & RM[l])==RM[l]))
			{
				genvalm_s[i] *= (1.0 + (s[k][l]*hs[k][l]));
				if (rnd == 1)	/*11*/;
				else
				{
					genvalm_s_F[i] *= (1.0 + s[k][l]);
				}
			}
		}

		for (k=0; k<NCRO; k++)
		for (l=1; l<NLOCI; l++)
		if ((initialgen[k][l] != (-99)) && (hs[k][l] == 99))
		{
			if (((gm[i][k][0] & RM[l])==RM[l])&&((gm[i][k][1] & RM[l])==RM[l]))	genvalm_s[i] *= (1.0 - sod2);
			else if (((gm[i][k][0] & RM[l])!=RM[l])&&((gm[i][k][1] & RM[l])!=RM[l]))	genvalm_s[i] *= (1.0 - sod1);
			else	
			{
				genvalm_s[i] *= 1.0;
				numAaSD ++;
			}
		}

	 	if (genvalm_s[i]<0.0)	genvalm_s[i]=0.0;

//		if (tracelevel!=0)	fprintf(fptr,"i = %d   addedfixed_s = %f  sod1 = %f  sod2 = %f  numAaSD = %f\n", i, addedfixed_s, sod1, sod2, numAaSD);
	}

	/* ******************* ADDITIVE AND DOMINANCE VARIANCE AND INBREEDING DEPRESSION RATE ***************** */

	va_del = 0.0;
	vd_del = 0.0;
	id_del = 0.0;
    
	va_let = 0.0;
	vd_let = 0.0;
	id_let = 0.0;
    
	va_od = 0.0;
	vd_od = 0.0;
	id_od = 0.0;
    
	numSegSD = 0;

	for (k=0; k<NCRO; k++)
	for (l=0; l<NLOCI; l++)
	if ((q[k][l] != 0.0)&&(q[k][l] != 1.0)&&(s[k][l] != 0.0))
	{
		if (s[k][l]<= (-0.9))
		{
			d_a_let = (s[k][l]/2.0) * (2.0*hs[k][l] - 1.0);
			alfa_a_let = (-s[k][l]/2.0) + ( d_a_let * (2.0*q[k][l] - 1.0) );
			va_let += 2.0 * q[k][l] * (1.0-q[k][l]) * alfa_a_let * alfa_a_let;
			vd_let += (2.0 * d_a_let * q[k][l] * (1.0-q[k][l])) * (2.0 * d_a_let * q[k][l] * (1.0-q[k][l]));
			id_let += (2.0 * d_a_let * q[k][l] * (1.0-q[k][l]));
		}
		else if (hs[k][l]==99)
		{
			numSegSD ++;

			d_a_od = (sod1 + sod2) / 2.0;
			alfa_a_od = (q[k][l]*sod2) - ((1.0-q[k][l])*sod1);
			va_od += 2.0 * alfa_a_od * alfa_a_od * q[k][l] * (1.0 - q[k][l]);
			vd_od += (2.0 * d_a_od * q[k][l] * (1.0-q[k][l])) * (2.0 * d_a_od * q[k][l] * (1.0-q[k][l]));
			id_od += (2.0 * d_a_od * q[k][l] * (1.0-q[k][l]));
			if (tracelevel!=0)   fprintf(fptr,"gen=%d k=%d l=%d q=%lf d=%lf sod1=%lf sod2=%lf va=%lf vd=%lf id=%lf\n", gen, k, l, q[k][l], d_a_od, sod1, sod2, va_od, vd_od, id_od);
		}
		else
		{
			d_a_del = (s[k][l]/2.0) * (2.0*hs[k][l] - 1.0);
			alfa_a_del = (-s[k][l]/2.0) + ( d_a_del * (2.0*q[k][l] - 1.0) );
 			va_del += 2.0 * q[k][l] * (1.0-q[k][l]) * alfa_a_del * alfa_a_del;
 			vd_del += (2.0 * d_a_del * q[k][l] * (1.0-q[k][l])) * (2.0 * d_a_del * q[k][l] * (1.0-q[k][l]));
			id_del += (2.0 * d_a_del * q[k][l] * (1.0-q[k][l]));
		}
	}

	accum (&AVE_VA_let[gen], va_let);
	accum (&AVE_VD_let[gen], vd_let);
	accum (&AVE_ID_let[gen], id_let);

	accum (&AVE_VA_od[gen], va_od);
	accum (&AVE_VD_od[gen], vd_od);
	accum (&AVE_ID_od[gen], id_od);

	accum (&AVE_VA_del[gen], va_del);
	accum (&AVE_VD_del[gen], vd_del);
	accum (&AVE_ID_del[gen], id_del);

	accum (&AVE_VA[gen], va_let + va_od + va_del);
	accum (&AVE_VD[gen], vd_let + vd_od + vd_del);
	accum (&AVE_ID[gen], id_let + id_od + id_del);
    
	accum (&Wmax[gen], addedfixed_s);
}


/* ***************************************************** */


void real_values (double s[][31], double hs[][31], double genvalm_s[])
{
	int AA, Aa, aa;
	double valor_AA, valor_Aa, valor_aa;
 
	// DELETEREOS
	for (k=0; k<NCRO; k++)
	for (l=1; l<NLOCI; l++)
	if ((initialgen[k][l] != (-99)) && (q[k][l] != 0.0) && (q[k][l] != 1.0) && (s[k][l] > (-0.9)) && (hs[k][l] != 99))
	{
		s_real[k][l] = 0.0;

		valor_AA=0.0; valor_Aa=0.0; valor_aa=0.0;
		AA = 0; Aa = 0; aa = 0;
		for (i=0; i<NIND; i++)
		{
			//if (tracelevel!=0)   fprintf(fptr,"\ni=%d genvalm_s=%f\n", i, genvalm_s[i]);

			if (((gm[i][k][0] & RM[l])==RM[l])&&((gm[i][k][1] & RM[l])==RM[l]))
			{
				aa ++;
				valor_aa += genvalm_s[i];
			}
			else if (((gm[i][k][0] & RM[l])!=RM[l])&&((gm[i][k][1] & RM[l])!=RM[l]))
			{
				AA ++;
				valor_AA += genvalm_s[i];
			}
			else
			{
				Aa ++;
				valor_Aa += genvalm_s[i];
			}
		}
		valor_AA = valor_AA / (double) AA;
		valor_Aa = valor_Aa / (double) Aa;
		valor_aa = valor_aa / (double) aa;

		if (tracelevel!=0)   if (aa != 0)	fprintf(fptr,"\nk=%d l=%d AA=%d vAA=%f Aa=%d vAa=%f aa=%d vaa=%f s=%f hs=%f q=%f\n",
						k, l, AA, valor_AA, Aa, valor_Aa, aa, valor_aa, s[k][l], hs[k][l], q[k][l]);
		valor_Aa = valor_Aa / valor_AA;
		valor_aa = valor_aa / valor_AA;
		valor_AA = 1.0;

		if ((AA != 0)&&(Aa != 0)&&(aa != 0))
		{
			s_real[k][l] = 1.0 - valor_aa;
			if (s_real[k][l] != 0.0) hs_real[k][l] = (1.0 - valor_Aa) / s_real[k][l];

			s_real[k][l] = -s_real[k][l];
			accum (&sdel_real[gen], s_real[k][l]);
			accum (&hsdel_real[gen], hs_real[k][l]);

			accum (&sdel_X[gen], s[k][l]);
			accum (&hsdel_X[gen], hs[k][l]);
		}

		if (tracelevel!=0)   if (aa != 0)	fprintf(fptr,"\nk=%d l=%d vsAA=%f vsAa=%f vsaa=%f s_real=%f hs_real=%f\n",
						k, l, valor_AA, valor_Aa, valor_aa, s_real[k][l], hs_real[k][l]);
	}


	/* ******************* ADDITIVE AND DOMINANCE VARIANCE AND INBREEDING DEPRESSION RATE OF EFFECTIVE VALUES ***************** */

	va_del = 0.0;
	vd_del = 0.0;
	id_del = 0.0;

	for (k=0; k<NCRO; k++)
	for (l=0; l<NLOCI; l++)
	if ((q[k][l] != 0.0)&&(q[k][l] != 1.0)&&(s_real[k][l] != 0.0)&&(hs[k][l] != 99))
	{
		d_a_del = (s_real[k][l]/2.0) * (2.0*hs_real[k][l] - 1.0);
		alfa_a_del = (-s_real[k][l]/2.0) + ( d_a_del * (2.0*q[k][l] - 1.0) );
		va_del += 2.0 * q[k][l] * (1.0-q[k][l]) * alfa_a_del * alfa_a_del;
		vd_del += (2.0 * d_a_del * q[k][l] * (1.0-q[k][l])) * (2.0 * d_a_del * q[k][l] * (1.0-q[k][l]));
		id_del += (2.0 * d_a_del * q[k][l] * (1.0-q[k][l]));
	}

	accum (&AVE_VA_del_real[gen], va_del);
	accum (&AVE_VD_del_real[gen], vd_del);
	accum (&AVE_ID_del_real[gen], id_del);

	id_del = 0.0;
	for (k=0; k<NCRO; k++)
	for (l=0; l<NLOCI; l++)
	if ((q[k][l] != 0.0)&&(q[k][l] != 1.0)&&(s_real[k][l] != 0.0)&&(hs[k][l] != 99))
	{
		d_a_del = (s[k][l]/2.0) * (2.0*hs[k][l] - 1.0);
		id_del += (2.0 * d_a_del * q[k][l] * (1.0-q[k][l]));
	}

	accum (&AVE_ID_del_X[gen], id_del);

//	SOBREDOMINANTES

	for (k=0; k<NCRO; k++)
	for (l=0; l<NLOCI; l++)
	if ((initialgen[k][l] != (-99)) && (q[k][l] != 0.0) && (q[k][l] != 1.0) && (hs[k][l] == 99))
	{
		s1_real[k][l] = 0.0;
		s2_real[k][l] = 0.0;

		valor_AA=0.0; valor_Aa=0.0; valor_aa=0.0;
		AA = 0; Aa = 0; aa = 0;
		for (i=0; i<NIND; i++)
		{
			//if (tracelevel!=0)   fprintf(fptr,"\ni=%d genvalm_s=%f\n", i, genvalm_s[i]);

			if (((gm[i][k][0] & RM[l])==RM[l])&&((gm[i][k][1] & RM[l])==RM[l]))
			{
				aa ++;
				valor_aa += genvalm_s[i];
			}
			else if (((gm[i][k][0] & RM[l])!=RM[l])&&((gm[i][k][1] & RM[l])!=RM[l]))
			{
				AA ++;
				valor_AA += genvalm_s[i];
			}
			else
			{
				Aa ++;
				valor_Aa += genvalm_s[i];
			}
		}
		valor_AA = valor_AA / (double) AA;
		valor_Aa = valor_Aa / (double) Aa;
		valor_aa = valor_aa / (double) aa;

		if (tracelevel!=0)   if (aa != 0)	fprintf(fptr,"\nk=%d l=%d AA=%d vAA=%f Aa=%d vAa=%f aa=%d vaa=%f s=%f hs=%f q=%f\n",
						k, l, AA, valor_AA, Aa, valor_Aa, aa, valor_aa, s[k][l], hs[k][l], q[k][l]);
		valor_AA = valor_AA / valor_Aa;
		valor_aa = valor_aa / valor_Aa;
		valor_Aa = 1.0;

		if ((AA != 0)&&(Aa != 0)&&(aa != 0))
		{
			s1_real[k][l] = 1.0 - valor_AA;
			s2_real[k][l] = 1.0 - valor_aa;

			accum (&sod1_real[gen], s1_real[k][l]);
			accum (&sod2_real[gen], s2_real[k][l]);
		}

		if (tracelevel!=0)   if (aa != 0)	fprintf(fptr,"\nk=%d l=%d vsAA=%f vsAa=%f vsaa=%f s1_real=%f s2_real=%f\n",
						k, l, valor_AA, valor_Aa, valor_aa, s1_real[k][l], s2_real[k][l]);
	}
	va_od = 0.0;
	vd_od = 0.0;
	id_od = 0.0;

	for (k=0; k<NCRO; k++)
	for (l=0; l<NLOCI; l++)
	if ((q[k][l] != 0.0)&&(q[k][l] != 1.0)&&(s1_real[k][l] != 0.0)&&(s2_real[k][l] != 0.0)&&(hs[k][l] == 99))
	{
		d_a_od = (s1_real[k][l] + s2_real[k][l]) / 2.0;
		alfa_a_od = (q[k][l]*s2_real[k][l]) - ((1.0-q[k][l])*s1_real[k][l]);
		va_od += 2.0 * alfa_a_od * alfa_a_od * q[k][l] * (1.0 - q[k][l]);
		vd_od += (2.0 * d_a_od * q[k][l] * (1.0-q[k][l])) * (2.0 * d_a_od * q[k][l] * (1.0-q[k][l]));
		id_od += (2.0 * d_a_od * q[k][l] * (1.0-q[k][l]));
		//if (tracelevel!=0) if (gen == 1)   fprintf(fptr,"gen=%d k=%d l=%d q=%lf d=%lf sod1=%f sod2=%f va=%f vd=%f id=%f\n", gen, k, l, q[k][l], d_a_od, sod1, sod2, va_od, vd_od, id_od);
	}
	accum (&AVE_VA_od_real[gen], va_od);
	accum (&AVE_VD_od_real[gen], vd_od);
	accum (&AVE_ID_od_real[gen], id_od);

	id_od = 0.0;
	for (k=0; k<NCRO; k++)
	for (l=0; l<NLOCI; l++)
	if ((q[k][l] != 0.0)&&(q[k][l] != 1.0)&&(s1_real[k][l] != 0.0)&&(s2_real[k][l] != 0.0)&&(hs[k][l] == 99))
	{
		d_a_od = (sod1 + sod2) / 2.0;
		id_od += (2.0 * d_a_od * q[k][l] * (1.0-q[k][l]));
	}
	accum (&AVE_ID_od_X[gen], id_od);
}


/* ***************************************************** */


void phenotypeB (double genvalm_s[], double  pm_s[])
{
	int ii, it;
	double gsum_s=0.0, gsum2_s=0.0;
    double gsum_s_F=0.0, gsum2_s_F=0.0;
	double maxpm=0.0;

	for (i=0; i<NIND; i++)
	{
		pm_s[i]=genvalm_s[i];

		gsum_s += pm_s[i];
		gsum2_s += (pm_s[i]*pm_s[i]);
        
        pm_s_F[i]=genvalm_s_F[i];

		gsum_s_F += pm_s_F[i];
		gsum2_s_F += (pm_s_F[i]*pm_s_F[i]);
	}

	accum (&gmean_s[gen], gsum_s/(double)NIND);
    accum (&gmean_s_F[gen], gsum_s_F/(double)NIND);
    if (gsum_s_F > 0.0)  accum (&ID_HOM[gen], -log(gsum_s_F/gsum_s));
	accum (&gvar_s[gen], (gsum2_s - (gsum_s*gsum_s / (double)NIND)) / ((double)NIND - 1.0));

	//if (tracelevel!=0)   fprintf(fptr,"\ngen = %d gmean_s = %f  gmean_s_F = %f ID_hom = %f\n", gen, gsum_s/(double)NIND, gsum_s_F/(double)NIND, -log(gsum_s_F/gsum_s));

	/* find individual with maximal fitness */
	//ii=0;
	//for (it=1; it<NIND; it++)	if (pm_s[ii] < pm_s[it])	ii=it;
	//maxpm = pm_s[ii];

	//for (i=0; i<NIND; i++)	    pm_s[i] = pm_s[i] / maxpm;
}


/* ***************************************************** */


dumpphenotypes()
{
	//if (tracelevel==0)   return (0);

	//fprintf(fptr,"\n Phenotypic value and Relative fitnesses \n");
	//for (i=0; i<NIND; i++)   fprintf(fptr,"%f\n", pm_s[i]);
}


/* ***************************************************** */

generation1 ()
{
	int a, p1, p2, EE[MM], FF[MM], family[NN], rnd;
	double family_sum=0.0, family_sum2=0.0;
	int numberrecs, nr, pointrec[MM][31], ncrorec[MM], rndk, rndl;

	for (i=0; i<NIND; i++)
	for (k=0; k<NCRO; k++)
	{
		sm[i][k][0]=gm[i][k][0];
		sm[i][k][1]=gm[i][k][1];
	}

	/***** disorder parents before mating *****/

	for (i=0; i<NIND; i++)
	{
		ran_i=(int)(uniform()*NIND);
		for (k=0; k<NCRO; k++)
		{
			a=sm[i][k][0]; sm[i][k][0]=sm[ran_i][k][0]; sm[ran_i][k][0]=a;
			a=sm[i][k][1]; sm[i][k][1]=sm[ran_i][k][1]; sm[ran_i][k][1]=a;
		}
	}

	/* *****Genotypic value of parents***** */

	//if (tracelevel!=0)	fprintf(fptr,"\n Parents of additional experiment \n");

	for (i=0; i<NIND; i++)
	{
		numAaSD = 0.0;

		pm_sG0[i] = addedfixed_s;

		// DELETEREOS
		for (k=0; k<NCRO; k++)
		for (l=1; l<NLOCI; l++)
		if ((initialgen[k][l] != (-99)) && (hs[k][l] != 99))
		{
			if (((sm[i][k][0] & RM[l])==RM[l])&&((sm[i][k][1] & RM[l])==RM[l]))
			{
				pm_sG0[i] *= (1.0 + s[k][l]);
			}
			else    if (((sm[i][k][0] & RM[l])!=RM[l])&&((sm[i][k][1] & RM[l])!=RM[l])) /* AA */;
			else
			{
			 	pm_sG0[i] *= (1.0 + (s[k][l]*hs[k][l]));
			}
		}

		for (k=0; k<NCRO; k++)
		for (l=1; l<NLOCI; l++)
		if ((initialgen[k][l] != (-99)) && (hs[k][l] == 99))
		{
			if (((sm[i][k][0] & RM[l])==RM[l])&&((sm[i][k][1] & RM[l])==RM[l]))		pm_sG0[i] *= (1.0 - sod2);
			else if (((sm[i][k][0] & RM[l])!=RM[l])&&((sm[i][k][1] & RM[l])!=RM[l]))	pm_sG0[i] *= (1.0 - sod1);
			else	
			{
                 pm_sG0[i] *= 1.0;
                 numAaSD ++;
			}
		}

		//if (tracelevel!=0)	fprintf(fptr," %d    sod1 = %f  sod2 = %f  numAaSD = %f  pm_sG0 = %f\n", i, sod1, sod2, numAaSD, pm_sG0[i]);
	}

	/******************************************/

	//if (tracelevel!=0)	fprintf(fptr,"\n Families of additional experiment \n");

	for (i=0; i<NIND;i++)   family[i]=0;

	for (i=0; i<((NIND/2)*J); i++)
	{
		/***** contributions from pairs *****/

		p1 = (i/J)*2;
		p2 = (i/J)*2 + 1;

		/*************************/

		//if (tracelevel!=0)   fprintf (fptr,"%d\t%d\n", p1, p2);

	    	if(L==99.0)
	    	{	    /* ******************* Free recombination ******************* */
			for (k=0; k<NCRO; k++)
			{
			   	EE[k] = (int)(uniform()*(pow(2.0,(double)NLOCI)));
				FF[k] = ~EE[k];
			   	gmG1[i][k][0]=((EE[k]&sm[p1][k][0])|(FF[k]&sm[p1][k][1]));
			}
			// if (tracelevel!=0)   fprintf (fptr,"i=%d EE[0]=%d EE[1]=%d EE[2]=%d sm00=%d sm01=%d sm10=%d sm11=%d sm20=%d sm21=%d g00=%d g10=%d g20=%d \n", i, 			EE[0], EE[1], EE[2], sm[p1][0][0], sm[p1][0][1], sm[p1][1][0], sm[p1][1][1], 				sm[p1][2][0], sm[p1][2][1], gmG1[i][0][0], gmG1[i][1][0], gmG1[i][2][0]);

			for (k=0; k<NCRO; k++)
			{
			   	EE[k] = (int)(uniform()*(pow(2.0,(double)NLOCI)));
			   	FF[k] = ~EE[k];
			   	gmG1[i][k][1]=((EE[k]&sm[p2][k][0])|(FF[k]&sm[p2][k][1]));
			}
			// if (tracelevel!=0)   fprintf (fptr,"i=%d EE[0]=%d EE[1]=%d EE[2]=%d sm00=%d sm01=%d sm10=%d sm11=%d sm20=%d sm21=%d g00=%d g10=%d g20=%d \n", i, 			EE[0], EE[1], EE[2], sm[p1][0][0], sm[p1][0][1], sm[p1][1][0], sm[p1][1][1], 				sm[p1][2][0], sm[p1][2][1], gmG1[i][0][0], gmG1[i][1][0], gmG1[i][2][0]);
		}
		else
		{	    /* ************** Restricted recombination ***************** */

			/* ****** Chromosome from father ****** */

			for (k=0; k<NCRO; k++)
			{
				ncrorec[k] = 0;
				for (l=0; l<NLOCI; l++)  pointrec[k][l] = 0;
			}

			// CROSSINGOVERS

			numberrecs = recombinationnumber();
			accum(&numberofrecom, numberrecs);
//			if (tracelevel!=0)   fprintf (fptr,"numberrecs=%d\n",numberrecs);

			for (nr=0; nr<numberrecs; nr++)
			{
				rndk = (int)(uniform()*NCRO);
				rndl = (int)(uniform()*NLOCI);
				ncrorec[rndk] = 1;
				pointrec[rndk][rndl] = 1;
//				if (tracelevel!=0)	fprintf (fptr,"Rec %d rndk=%d  rndl=%d\n", nr, rndk, rndl);
			}

			marker = 1;

			for (k=0; k<NCRO; k++)
			{
				EE[k]=0;
				if (ncrorec[k] == 0)
				{
					if (marker==(-1))
					{
						EE[k] = ~EE[k];
					}
				}
				else
				{
					for (l=0; l<NLOCI; l++)
			      	{
						if (pointrec[k][l] == 0)
						{
							if (marker==(-1))  EE[k] = EE[k] | RM[l];
						}
						else
						{
							if (marker==1)
							{
								EE[k] = EE[k] | RM[l];
								marker = marker * (-1);
							}
							else
							{
								marker = marker * (-1);
							}
						}
					}
				}
			}

			rnd = uniform();
			for (k=0; k<NCRO; k++)
			{
				if (rnd < 0.5)	EE[k] = ~EE[k];
				FF[k] = ~EE[k];
				gmG1[i][k][0]=((EE[k]&sm[p1][k][0])|(FF[k]&sm[p1][k][1]));
			}

			/* ****** Chromosome from mother ****** */

			for (k=0; k<NCRO; k++)
			{
				ncrorec[k] = 0;
				for (l=0; l<NLOCI; l++)  pointrec[k][l] = 0;
			}

			// CROSSINGOVERS

			numberrecs = recombinationnumber();
			accum(&numberofrecom, numberrecs);
//			if (tracelevel!=0)   fprintf (fptr,"numberrecs=%d\n",numberrecs);

			for (nr=0; nr<numberrecs; nr++)
			{
				rndk = (int)(uniform()*NCRO);
				rndl = (int)(uniform()*NLOCI);
				ncrorec[rndk] = 1;
				pointrec[rndk][rndl] = 1;
//				if (tracelevel!=0)	fprintf (fptr,"Rec %d rndk=%d  rndl=%d\n", nr, rndk, rndl);
			}

			marker = 1;

			for (k=0; k<NCRO; k++)
			{
				EE[k]=0;
				if (ncrorec[k] == 0)
				{
					if (marker==(-1))
					{
						EE[k] = ~EE[k];
					}
				}
				else
				{
					for (l=0; l<NLOCI; l++)
			      	{
						if (pointrec[k][l] == 0)
						{
							if (marker==(-1))  EE[k] = EE[k] | RM[l];
						}
						else
						{
							if (marker==1)
							{
								EE[k] = EE[k] | RM[l];
								marker = marker * (-1);
							}
							else
							{
								marker = marker * (-1);
							}
						}
					}
				}
			}

			rnd = uniform();
			for (k=0; k<NCRO; k++)
			{
				if (rnd < 0.5)	EE[k] = ~EE[k];
				FF[k] = ~EE[k];
				gmG1[i][k][1]=((EE[k]&sm[p2][k][0])|(FF[k]&sm[p2][k][1]));

//				if (tracelevel!=0)	if (k<=3) if ((EE[k] & RM[0]) == RM[0]) fprintf(fptr,"EE[%d]=1\n", k);
			}
		}

		/* *******************PPPPPPPPPPPPPP******************** */

		/* *****Genotypic value of offspring***** */

		numAaSD = 0.0;

		pm_sG1[i] = addedfixed_s;

		// DELETEREOS
		for (k=0; k<NCRO; k++)
		for (l=1; l<NLOCI; l++)
		if ((initialgen[k][l] != (-99)) && (hs[k][l] != 99))
		{
			if (((gmG1[i][k][0] & RM[l])==RM[l])&&((gmG1[i][k][1] & RM[l])==RM[l]))
			{
				pm_sG1[i] *= (1.0 + s[k][l]);
			}
			else    if (((gmG1[i][k][0] & RM[l])!=RM[l])&&((gmG1[i][k][1] & RM[l])!=RM[l])) /* AA */;
			else
			{
			 	pm_sG1[i] *= (1.0 + (s[k][l]*hs[k][l]));
			}
		}

		for (k=0; k<NCRO; k++)
		for (l=1; l<NLOCI; l++)
		if ((initialgen[k][l] != (-99)) && (hs[k][l] == 99))
		{
			if (((gmG1[i][k][0] & RM[l])==RM[l])&&((gmG1[i][k][1] & RM[l])==RM[l]))		pm_sG1[i] *= (1.0 - sod2);
			else if (((gmG1[i][k][0] & RM[l])!=RM[l])&&((gmG1[i][k][1] & RM[l])!=RM[l]))	pm_sG1[i] *= (1.0 - sod1);
			else	
			{
                 pm_sG1[i] *= 1.0;
                 numAaSD ++;
			}
		}

		//if (tracelevel!=0)	fprintf(fptr," %d    sod1 = %f  sod2 = %f  numAaSD = %f  pm_sG1 = %f\n", i, sod1, sod2, numAaSD, pm_sG1[i]);
        
		family[p1]+=1;    family[p2]+=1;
	}

	/* ***** VARIANCE OF FAMILY SIZE ***** */

	for (i=0; i<NIND; i++)
	{
		//if (tracelevel!=0)	fprintf(fptr," family[%d] = %d\n", i, family[i]);
		family_sum += family[i];
		family_sum2 += (family[i] * family[i]);
	}

	//if (tracelevel!=0)	fprintf (fptr,"SK2 = %5.3f\n", ((family_sum2-(family_sum*family_sum/NIND))/(NIND-1.0)));

	/* ***** ESTIMATION OF ADDITIVE VARIANCE FROM COVARIANCE OF MEAN_OFFSPRING AND MEAN_PARENT ***** */

	for (i=0; i<NIND; i+=2)
	{
		mean_parents_G0[i] = 0.0;
		mean_family_G1[i] = 0.0;
	}

//	LOG SCALING
	I = NIND/2;
	for (i=0; i<NIND; i+=2)	if ((pm_sG0[i] != 0.0)&&(pm_sG0[i+1] != 0.0)) mean_parents_G0[i] += (log(pm_sG0[i]) + log(pm_sG0[i+1])) / 2.0;
	for (i=0; i<(I*J); i++)	if (pm_sG1[i] != 0.0) mean_family_G1[(i/J)*2] += log(pm_sG1[i]) / J;

	//if (tracelevel!=0)	for (i=0; i<NIND; i+=2)	fprintf (fptr,"mean_parents_G0[%d] = %f    mean_family_G1[%d] = %f\n", i, mean_parents_G0[i], i, mean_family_G1[i]);

	double sum_parents = 0.0;
	double sum_offspring = 0.0;
	double sum_parents_offspring = 0.0;
	double VA_G0;

	for (i=0; i<NIND; i+=2)
	{
		sum_parents += mean_parents_G0[i];
		sum_offspring += mean_family_G1[i];
		sum_parents_offspring += mean_parents_G0[i] * mean_family_G1[i];
	}
	VA_G0 = (sum_parents_offspring - (sum_parents * sum_offspring)/I) / (I - 1.0);
    
	//if (tracelevel!=0)	fprintf (fptr,"sum_parents = %f   sum_offspring = %f   sum_parents_offspring = %f\n", sum_parents, sum_offspring, sum_parents_offspring);
	//if (tracelevel!=0)	fprintf (fptr,"VA_G0 = %f\n", 2.0*VA_G0);

	accum(&AVE_VAeff[gen], 2.0*VA_G0);

	/* ***** ESTIMATION OF ADDITIVE AND DOMINANCE VARIANCE FROM ANOVA FULL-SIB FAMILIES ***** */

	double Xoo = 0.0, Xio2 = 0.0, Xij2 = 0.0;
	double SSF, SSE, MSF, MSE, vf;

//	LOG SCALING
	for (i=0; i<((NIND/2)*J); i++)
	{
		if (pm_sG1[i] != 0.0)	Xij2 += (log(pm_sG1[i]) * log(pm_sG1[i]));
		if (pm_sG1[i] != 0.0)	Xoo += log(pm_sG1[i]);
	}
	for (i=0; i<NIND; i+=2)		Xio2 += (J * mean_family_G1[i] * J * mean_family_G1[i]);
    
	//if (tracelevel!=0)	fprintf (fptr,"Xij2 = %f  Xoo = %f  Xio2 = %f\n", Xij2, Xoo, Xio2);


	SSF = (Xio2 / J) - (Xoo * Xoo / (I*J));
	SSE = Xij2 - (Xio2 / J);
    
	//if (tracelevel!=0)	fprintf (fptr,"SSF = %f  SSE = %f\n", SSF, SSE);

	MSF = SSF / (I - 1.0);
	MSE = SSE / (I*(J - 1.0));

	vf = (MSF - MSE) / J;

	//if (tracelevel!=0)	fprintf (fptr,"MSF = %f  MSE = %f  vf = %f    VD = %f\n", MSF, MSE, vf, 4.0*(vf - VA_G0));

	accum(&AVE_VDeff[gen], 4.0*(vf - VA_G0));

	/* ***** FREQUENCY OF NEUTRAL ALLELES IN PROGENY ***** */

	double qG1[MM][31];

	for (k=0; k<NCRO; k++)
	if(initialgen[k][0] != (-99))
	{
		AA=0.0; Aa=0.0; aa=0.0;

		for (i=0; i<(I*J); i++)
		{
			if (((gmG1[i][k][0] & RM[0])==RM[0])&&((gmG1[i][k][1] & RM[0])==RM[0]))		aa+=1.0;
			else if (((gmG1[i][k][0] & RM[0])!=RM[0])&&((gmG1[i][k][1] & RM[0])!=RM[0]))	AA+=1.0;
			else	Aa+=1.0;
		}

		qG1[k][0] = (aa/(double)(I*J))+(Aa/(2.0*(double)(I*J)));
	}

	/* ***** Estimates of F_YANG ***** */

	double numLOCI, sum_YANG2[NN];

	numLOCI = 0.0;

	for (k=0; k<NCRO; k++)
	if((qG1[k][0] != 0.0)&&(qG1[k][0] != 1.0))	numLOCI ++;

	for (i=0; i<(I*J); i++)
	{
		sum_YANG2[i] = 0.0; 

		for (k=0; k<NCRO; k++)
		if((qG1[k][0] != 0.0)&&(qG1[k][0] != 1.0))
		{
			if (((gmG1[i][k][0] & RM[0])==RM[0])&&((gmG1[i][k][1] & RM[0])==RM[0]))
				sum_YANG2[i] += ( (4.0)-((1.0+2.0*qG1[k][0])*2.0)+(2.0*qG1[k][0]*qG1[k][0]) ) / (2.0*qG1[k][0]*(1.0-qG1[k][0])); 
			else if (((gmG1[i][k][0] & RM[0])!=RM[0])&&((gmG1[i][k][1] & RM[0])!=RM[0]))
				sum_YANG2[i] += (2.0*qG1[k][0]*qG1[k][0]) / (2.0*qG1[k][0]*(1.0-qG1[k][0])); 
			else	sum_YANG2[i] += ( (1.0)-((1.0+2.0*qG1[k][0])*1.0)+(2.0*qG1[k][0]*qG1[k][0]) ) / (2.0*qG1[k][0]*(1.0-qG1[k][0])); 
		}

		fyang2[i] = sum_YANG2[i] / numLOCI;

//		if (tracelevel!=0) fprintf(fptr,"gen=%d ind=%d numLOCI=%f fyang2=%f\n", gen, i, numLOCI, fyang2[i]);
	}

	/* ***** INBREEDING DEPRESSION FROM F_YANG ***** */

	double sum_G, sum_F, sum_F2, sum_FG, var_F;

	sum_G=0.0; sum_F=0.0; sum_F2=0.0; sum_FG=0.0;
	for (i=0; i<(I*J); i++)
	{
		if (pm_sG1[i] == 0.0) pm_sG1[i]=0.000001;

		sum_F += fyang2[i];
		sum_F2 += (fyang2[i] * fyang2[i]);
		sum_G += log(pm_sG1[i]);
		sum_FG += (fyang2[i] * log(pm_sG1[i]));

		//if (tracelevel!=0) fprintf(fptr,"gen=%d ind=%d pm_sG1=%f fyang2=%f\n", gen, i, pm_sG1[i], fyang2[i]);
	}
	var_F = sum_F2 - (sum_F * sum_F / (double)(I*J));
	if (var_F > 0.0) accum (&IDFyang2[gen], (sum_FG - (sum_G * sum_F / (double)(I*J))) / var_F);

	//if (tracelevel!=0) if (var_F > 0.0) fprintf(fptr,"*****gen=%d ID=%f\n", gen, (sum_FG - (sum_G * sum_F / (double)(I*J))) / var_F);
}


/* ***************************************************** */


generation2 ()
{
	int a, p1, p2, EE[MM], FF[MM], family[NN], rnd;
	double family_sum=0.0, family_sum2=0.0;
	int numberrecs, nr, pointrec[MM][31], ncrorec[MM], rndk, rndl;

	/***** choose as parents first and second offspring from each family *****/

	//if (tracelevel!=0)	fprintf(fptr,"\n FS pairs of parents \n");

	for (i=0; i<NIND; i+=2)
	for (k=0; k<NCRO; k++)
	{
		sm[i][k][0]=gmG1[(i*J/2)][k][0];
		sm[i][k][1]=gmG1[(i*J/2)][k][1];

		sm[i+1][k][0]=gmG1[(i*J/2)+1][k][0];
		sm[i+1][k][1]=gmG1[(i*J/2)+1][k][1];

		//if ((tracelevel!=0)&&(k==0))	fprintf(fptr,"\n parent %d from Full sibs pm_sG1[%d]=%f and pm_sG1[%d]=%f\n", i, i*J/2, pm_sG1[i*J/2], (i*J/2)+1, pm_sG1[(i*J/2)+1]);
	}

	/******************************************/

	//if (tracelevel!=0)	fprintf(fptr,"\n Families from FS parents of additional experiment \n");

	for (i=0; i<NIND;i++)   family[i]=0;

	for (i=0; i<((NIND/2)*J); i++)
	{
		/***** contributions from pairs *****/

		p1 = (i/J)*2;
		p2 = (i/J)*2 + 1;

		/*************************/

		//if (tracelevel!=0)   fprintf (fptr,"%d\t%d\n", p1, p2);

	    	if(L==99.0)
	    	{	    /* ******************* Free recombination ******************* */
			for (k=0; k<NCRO; k++)
			{
			   	EE[k] = (int)(uniform()*(pow(2.0,(double)NLOCI)));
				FF[k] = ~EE[k];
			   	gmG2[i][k][0]=((EE[k]&sm[p1][k][0])|(FF[k]&sm[p1][k][1]));
			}
			// if (tracelevel!=0)   fprintf (fptr,"i=%d EE[0]=%d EE[1]=%d EE[2]=%d sm00=%d sm01=%d sm10=%d sm11=%d sm20=%d sm21=%d g00=%d g10=%d g20=%d \n", i, 			EE[0], EE[1], EE[2], sm[p1][0][0], sm[p1][0][1], sm[p1][1][0], sm[p1][1][1], 				sm[p1][2][0], sm[p1][2][1], gmG2[i][0][0], gmG2[i][1][0], gmG2[i][2][0]);

			for (k=0; k<NCRO; k++)
			{
			   	EE[k] = (int)(uniform()*(pow(2.0,(double)NLOCI)));
			   	FF[k] = ~EE[k];
			   	gmG2[i][k][1]=((EE[k]&sm[p2][k][0])|(FF[k]&sm[p2][k][1]));
			}
			// if (tracelevel!=0)   fprintf (fptr,"i=%d EE[0]=%d EE[1]=%d EE[2]=%d sm00=%d sm01=%d sm10=%d sm11=%d sm20=%d sm21=%d g00=%d g10=%d g20=%d \n", i, 			EE[0], EE[1], EE[2], sm[p1][0][0], sm[p1][0][1], sm[p1][1][0], sm[p1][1][1], 				sm[p1][2][0], sm[p1][2][1], gmG2[i][0][0], gmG2[i][1][0], gmG2[i][2][0]);
		}
		else
		{	    /* ************** Restricted recombination ***************** */

			/* ****** Chromosome from father ****** */

			for (k=0; k<NCRO; k++)
			{
				ncrorec[k] = 0;
				for (l=0; l<NLOCI; l++)  pointrec[k][l] = 0;
			}

			// CROSSINGOVERS

			numberrecs = recombinationnumber();
			accum(&numberofrecom, numberrecs);
//			if (tracelevel!=0)   fprintf (fptr,"numberrecs=%d\n",numberrecs);

			for (nr=0; nr<numberrecs; nr++)
			{
				rndk = (int)(uniform()*NCRO);
				rndl = (int)(uniform()*NLOCI);
				ncrorec[rndk] = 1;
				pointrec[rndk][rndl] = 1;
//				if (tracelevel!=0)	fprintf (fptr,"Rec %d rndk=%d  rndl=%d\n", nr, rndk, rndl);
			}

			marker = 1;

			for (k=0; k<NCRO; k++)
			{
				EE[k]=0;
				if (ncrorec[k] == 0)
				{
					if (marker==(-1))
					{
						EE[k] = ~EE[k];
					}
				}
				else
				{
					for (l=0; l<NLOCI; l++)
			      	{
						if (pointrec[k][l] == 0)
						{
							if (marker==(-1))  EE[k] = EE[k] | RM[l];
						}
						else
						{
							if (marker==1)
							{
								EE[k] = EE[k] | RM[l];
								marker = marker * (-1);
							}
							else
							{
								marker = marker * (-1);
							}
						}
					}
				}
			}

			rnd = uniform();
			for (k=0; k<NCRO; k++)
			{
				if (rnd < 0.5)	EE[k] = ~EE[k];
				FF[k] = ~EE[k];
				gmG2[i][k][0]=((EE[k]&sm[p1][k][0])|(FF[k]&sm[p1][k][1]));
			}

			/* ****** Chromosome from mother ****** */

			for (k=0; k<NCRO; k++)
			{
				ncrorec[k] = 0;
				for (l=0; l<NLOCI; l++)  pointrec[k][l] = 0;
			}

			// CROSSINGOVERS

			numberrecs = recombinationnumber();
			accum(&numberofrecom, numberrecs);
//			if (tracelevel!=0)   fprintf (fptr,"numberrecs=%d\n",numberrecs);

			for (nr=0; nr<numberrecs; nr++)
			{
				rndk = (int)(uniform()*NCRO);
				rndl = (int)(uniform()*NLOCI);
				ncrorec[rndk] = 1;
				pointrec[rndk][rndl] = 1;
//				if (tracelevel!=0)	fprintf (fptr,"Rec %d rndk=%d  rndl=%d\n", nr, rndk, rndl);
			}

			marker = 1;

			for (k=0; k<NCRO; k++)
			{
				EE[k]=0;
				if (ncrorec[k] == 0)
				{
					if (marker==(-1))
					{
						EE[k] = ~EE[k];
					}
				}
				else
				{
					for (l=0; l<NLOCI; l++)
			      	{
						if (pointrec[k][l] == 0)
						{
							if (marker==(-1))  EE[k] = EE[k] | RM[l];
						}
						else
						{
							if (marker==1)
							{
								EE[k] = EE[k] | RM[l];
								marker = marker * (-1);
							}
							else
							{
								marker = marker * (-1);
							}
						}
					}
				}
			}

			rnd = uniform();
			for (k=0; k<NCRO; k++)
			{
				if (rnd < 0.5)	EE[k] = ~EE[k];
				FF[k] = ~EE[k];
				gmG2[i][k][1]=((EE[k]&sm[p2][k][0])|(FF[k]&sm[p2][k][1]));

//				if (tracelevel!=0)	if (k<=3) if ((EE[k] & RM[0]) == RM[0]) fprintf(fptr,"EE[%d]=1\n", k);
			}
		}

		/* *******************PPPPPPPPPPPPPP******************** */

		/* *****Genotypic value of offspring from FS pairs***** */

		numAaSD = 0.0;

		pm_sG2[i] = addedfixed_s;

		// DELETEREOS
		for (k=0; k<NCRO; k++)
		for (l=1; l<NLOCI; l++)
		if ((initialgen[k][l] != (-99)) && (hs[k][l] != 99))
		{
			if (((gmG2[i][k][0] & RM[l])==RM[l])&&((gmG2[i][k][1] & RM[l])==RM[l]))
			{
				pm_sG2[i] *= (1.0 + s[k][l]);
			}
			else    if (((gmG2[i][k][0] & RM[l])!=RM[l])&&((gmG2[i][k][1] & RM[l])!=RM[l])) /* AA */;
			else
			{
			 	pm_sG2[i] *= (1.0 + (s[k][l]*hs[k][l]));
			}
		}

		for (k=0; k<NCRO; k++)
		for (l=1; l<NLOCI; l++)
		if ((initialgen[k][l] != (-99)) && (hs[k][l] == 99))
		{
			if (((gmG2[i][k][0] & RM[l])==RM[l])&&((gmG2[i][k][1] & RM[l])==RM[l]))		pm_sG2[i] *= (1.0 - sod2);
			else if (((gmG2[i][k][0] & RM[l])!=RM[l])&&((gmG2[i][k][1] & RM[l])!=RM[l]))	pm_sG2[i] *= (1.0 - sod1);
			else	
			{
                 pm_sG2[i] *= 1.0;
                 numAaSD ++;
			}
		}

//		if (tracelevel!=0)	fprintf(fptr," %d    sod1 = %f  sod2 = %f  numAaSD = %f  pm_sG2 = %f\n", i, sod1, sod2, numAaSD, pm_sG2[i]);
        
		family[p1]+=1;    family[p2]+=1;
	}

	/* ***** VARIANCE OF FAMILY SIZE ***** */

	for (i=0; i<NIND; i++)
	{
//		if (tracelevel!=0)	fprintf(fptr," family[%d] = %d\n", i, family[i]);
		family_sum += family[i];
		family_sum2 += (family[i] * family[i]);
	}

//	if (tracelevel!=0)	fprintf (fptr,"SK2 = %5.3f\n", ((family_sum2-(family_sum*family_sum/NIND))/(NIND-1.0)));

	/* ***** ESTIMATION OF INBREEDING LOAD (B) FROM FULL SIB FAMILIES ***** */

	double mean_Wo=0.0, mean_Wfs=0.0, B_FS;

	for (i=0; i<NIND; i+=2)
	{
		mean_family_G1[i] = 0.0;
		mean_family_G2[i] = 0.0;
	}

	for (i=0; i<(I*J); i++)	mean_family_G1[(i/J)*2] += pm_sG1[i] / J;
	for (i=0; i<(I*J); i++)	mean_family_G2[(i/J)*2] += pm_sG2[i] / J;

	for (i=0; i<NIND; i+=2)
	{
		mean_Wo += mean_family_G1[i] / I;
		mean_Wfs += mean_family_G2[i] / I;
	}

	//if (tracelevel!=0)	for (i=0; i<NIND; i+=2)	fprintf (fptr,"mean_family_G1[%d] = %f   mean_family_G2[%d] = %f\n", i, mean_family_G1[i], i, mean_family_G2[i]);

	B_FS = log(mean_Wo / mean_Wfs) / 0.25;

	//if (tracelevel!=0)	fprintf (fptr,"mean_Wo = %f   mean_Wfs = %f   B_FS = %f\n", mean_Wo, mean_Wfs, B_FS);

	accum(&AVE_B_FS[gen], B_FS);
}


/* ***************************************************** */

void mating (int gm[][MM][2], int sm[][MM][2], double  pm_s[])
{
	int a, p1, p2, EE[MM], FF[MM], family[NN], rnd;
	double family_sum2=0.0, cum=0.0, cum_pm_s[NN];
	int numberrecs, nr, pointrec[MM][31], ncrorec[MM], rndk, rndl;

	for (i=0; i<NIND; i++)
	for (k=0; k<NCRO; k++)
	{
		sm[i][k][0]=gm[i][k][0];
		sm[i][k][1]=gm[i][k][1];
	}

	/********* Accumulate ps values ***********/

	for (i=0; i<NIND;i++)   family[i]=0;

	for (i=0; i<NIND;i++)	if (pm_s[i] == 0.0)	pm_s[i] = 0.01;

	for (i=0; i<NIND; i++)
	{
		cum += pm_s[i];
		cum_pm_s[i] = cum;
	}
	for (i=0; i<NIND; i++)	cum_pm_s[i] = cum_pm_s[i] / cum;

	//if (tracelevel!=0)	for (i=0; i<NIND; i++)	fprintf (fptr," %d    pm_s = %f   cum_pm_s = %f\n", i, pm_s[i], cum_pm_s[i]);

	/******************************************/

	//if (tracelevel!=0)	fprintf(fptr,"\n Parents \n");

	for (i=0; i<NIND; i++)
	{
		/***** parents *****/

	    p1 = binarysearchcumulative(cum_pm_s, NIND);

	    if (PS==99)
	    {
			do { p2 = binarysearchcumulative(cum_pm_s, NIND); }
			while (p2 == p1);
	    }
	    else
	    {
			if (uniform() < PS)		p2=p1;
			else
			{
				do { p2 = binarysearchcumulative(cum_pm_s, NIND); }
				while (p2 == p1);
			}
	    }

	    family[p1]+=1;    family[p2]+=1;

		/*************************/

		//if (tracelevel!=0)   fprintf (fptr,"%d\t%d\n", p1, p2);

	    if(L==99.0)
	    {	    /* ******************* Free recombination ******************* */
			for (k=0; k<NCRO; k++)
			{
			   	EE[k] = (int)(uniform()*(pow(2.0,(double)NLOCI)));
				FF[k] = ~EE[k];
			   	gm[i][k][0]=((EE[k]&sm[p1][k][0])|(FF[k]&sm[p1][k][1]));
			}
			// if (tracelevel!=0)   fprintf (fptr,"i=%d EE[0]=%d EE[1]=%d EE[2]=%d sm00=%d sm01=%d sm10=%d sm11=%d sm20=%d sm21=%d g00=%d g10=%d g20=%d \n", i, 			EE[0], EE[1], EE[2], sm[p1][0][0], sm[p1][0][1], sm[p1][1][0], sm[p1][1][1], 				sm[p1][2][0], sm[p1][2][1], gm[i][0][0], gm[i][1][0], gm[i][2][0]);

			for (k=0; k<NCRO; k++)
			{
			   	EE[k] = (int)(uniform()*(pow(2.0,(double)NLOCI)));
			   	FF[k] = ~EE[k];
			   	gm[i][k][1]=((EE[k]&sm[p2][k][0])|(FF[k]&sm[p2][k][1]));
			}
			// if (tracelevel!=0)   fprintf (fptr,"i=%d EE[0]=%d EE[1]=%d EE[2]=%d sm00=%d sm01=%d sm10=%d sm11=%d sm20=%d sm21=%d g00=%d g10=%d g20=%d \n", i, 			EE[0], EE[1], EE[2], sm[p1][0][0], sm[p1][0][1], sm[p1][1][0], sm[p1][1][1], 				sm[p1][2][0], sm[p1][2][1], gm[i][0][0], gm[i][1][0], gm[i][2][0]);
		}
		else
		{	    /* ************** Restricted recombination ***************** */

			/* ****** Chromosome from father ****** */

			for (k=0; k<NCRO; k++)
			{
				ncrorec[k] = 0;
				for (l=0; l<NLOCI; l++)  pointrec[k][l] = 0;
			}

			// CROSSINGOVERS

			numberrecs = recombinationnumber();
			accum(&numberofrecom, numberrecs);
//			if (tracelevel!=0)   fprintf (fptr,"numberrecs=%d\n",numberrecs);

			for (nr=0; nr<numberrecs; nr++)
			{
				rndk = (int)(uniform()*NCRO);
				rndl = (int)(uniform()*NLOCI);
				ncrorec[rndk] = 1;
				pointrec[rndk][rndl] = 1;
//				if (tracelevel!=0)	fprintf (fptr,"Rec %d rndk=%d  rndl=%d\n", nr, rndk, rndl);
			}

/*			if (i == 0)
			for (k=0; k<NCRO; k++)
			for (l=0; l<NLOCI; l++)
			{
				if (pointrec[k][l] == 1)
				if (tracelevel!=0)	fprintf (fptr,"pointrec k=%d l=%d\n", k, l);
			}
*/
			marker = 1;

			for (k=0; k<NCRO; k++)
			{
				EE[k]=0;
				if (ncrorec[k] == 0)
				{
					if (marker==(-1))
					{
						EE[k] = ~EE[k];
					}
				}
				else
				{
					for (l=0; l<NLOCI; l++)
			      	{
						if (pointrec[k][l] == 0)
						{
							if (marker==(-1))  EE[k] = EE[k] | RM[l];
						}
						else
						{
							if (marker==1)
							{
								EE[k] = EE[k] | RM[l];
								marker = marker * (-1);
							}
							else
							{
								marker = marker * (-1);
							}
						}
					}
				}
			}

			rnd = uniform();
			for (k=0; k<NCRO; k++)
			{
				if (rnd < 0.5)	EE[k] = ~EE[k];
				FF[k] = ~EE[k];
				gm[i][k][0]=((EE[k]&sm[p1][k][0])|(FF[k]&sm[p1][k][1]));
			}
/*
			if (tracelevel!=0)
			if (i == 0)
			{
				fprintf (fptr,"EE\n");
				for (k=0; k<NCRO; k++)
				{
					for (l=0; l<NLOCI; l++)
					{
						if((EE[k]&RM[l])==RM[l])  	fprintf (fptr,"1 ");
						else			   	fprintf (fptr,"0 ");
					}
					fprintf (fptr,"\n");
				}
				fprintf (fptr,"FF\n");
				for (k=0; k<NCRO; k++)
				{
					for (l=0; l<NLOCI; l++)
					{
						if((FF[k]&RM[l])==RM[l])  	fprintf (fptr,"1 ");
						else			   	fprintf (fptr,"0 ");
					}
					fprintf (fptr,"\n");
				}
			}
*/
			/* ****** Chromosome from mother ****** */

			for (k=0; k<NCRO; k++)
			{
				ncrorec[k] = 0;
				for (l=0; l<NLOCI; l++)  pointrec[k][l] = 0;
			}

			// CROSSINGOVERS

			numberrecs = recombinationnumber();
			accum(&numberofrecom, numberrecs);
//			if (tracelevel!=0)   fprintf (fptr,"numberrecs=%d\n",numberrecs);

			for (nr=0; nr<numberrecs; nr++)
			{
				rndk = (int)(uniform()*NCRO);
				rndl = (int)(uniform()*NLOCI);
				ncrorec[rndk] = 1;
				pointrec[rndk][rndl] = 1;
//				if (tracelevel!=0)	fprintf (fptr,"Rec %d rndk=%d  rndl=%d\n", nr, rndk, rndl);
			}

			marker = 1;

			for (k=0; k<NCRO; k++)
			{
				EE[k]=0;
				if (ncrorec[k] == 0)
				{
					if (marker==(-1))
					{
						EE[k] = ~EE[k];
					}
				}
				else
				{
					for (l=0; l<NLOCI; l++)
			      	{
						if (pointrec[k][l] == 0)
						{
							if (marker==(-1))  EE[k] = EE[k] | RM[l];
						}
						else
						{
							if (marker==1)
							{
								EE[k] = EE[k] | RM[l];
								marker = marker * (-1);
							}
							else
							{
								marker = marker * (-1);
							}
						}
					}
				}
			}

			rnd = uniform();
			for (k=0; k<NCRO; k++)
			{
				if (rnd < 0.5)	EE[k] = ~EE[k];
				FF[k] = ~EE[k];
				gm[i][k][1]=((EE[k]&sm[p2][k][0])|(FF[k]&sm[p2][k][1]));

//				if (tracelevel!=0)	if (k<=3) if ((EE[k] & RM[0]) == RM[0]) fprintf(fptr,"EE[%d]=1\n", k);
			}
		}

		/* *******************PPPPPPPPPPPPPP******************** */

		/* *****Genotypic value of offspring***** */

		pm_s[i]=addedfixed_s;

		// DELETEREOS
		for (k=0; k<NCRO; k++)
		for (l=1; l<NLOCI; l++)
		if ((initialgen[k][l] != (-99)) && (hs[k][l] != 99))
		{
			if (((gm[i][k][0] & RM[l])==RM[l])&&((gm[i][k][1] & RM[l])==RM[l]))
			{
				pm_s[i] *= (1.0 + s[k][l]);
			}
			else    if (((gm[i][k][0] & RM[l])!=RM[l])&&((gm[i][k][1] & RM[l])!=RM[l])) /* AA */;
			else
			{
			 	pm_s[i] *= (1.0 + (s[k][l]*hs[k][l]));
			}
		}

        //SD
        
        for (k=0; k<NCRO; k++)
		for (l=1; l<NLOCI; l++)
		if ((initialgen[k][l] != (-99)) && (hs[k][l] == 99))
		{
			if (((gm[i][k][0] & RM[l])==RM[l])&&((gm[i][k][1] & RM[l])==RM[l]))			pm_s[i] *= (1.0 - sod2);
			else if (((gm[i][k][0] & RM[l])!=RM[l])&&((gm[i][k][1] & RM[l])!=RM[l]))	pm_s[i] *= (1.0 - sod1);
			else		pm_s[i] *= 1.0;
		}
        

		//if (tracelevel!=0)	fprintf (fptr," %d    pm_s = %f\n", i, pm_s[i]);
	}

	/* ***** VARIANCE OF FAMILY SIZE ***** */

	for (i=0; i<NIND; i++)	family_sum2 += (family[i] * family[i]);

	accum(&SK2[gen], ((family_sum2-(4.0*NIND))/(NIND-1.0)) );
	//if (tracelevel!=0)	fprintf (fptr,"SK2 = %5.3f\n", ((family_sum2-(4.0*NIND))/(NIND-1.0)));
}


/* ***************************************************** */


int recombinationnumber ()
{
	int r;
	if ((L < normalthreshold) && (exp(-L) != 0.0) )
	{
		r = poisson(lastinrecpoissontable, recpoissontable);
	}
	else r = (int)normal(L, sqrt(L));
	return(r);
}


/* ***************************************************** */


dumpoffspring()
{
	//if (tracelevel==0)   return (0);

	//fprintf(fptr,"\n Offspring before mutation (gm0 gm1)\n");
	//for (i=0; i<NIND; i++)   fprintf(fptr,"%d  %d\n",gm[i][0][0],gm[i][0][1]);
}


/* ***************************************************** */


printout()
{
	fgen = fopen ("genfile.dat","a");
	for (gen=0; gen<=generations; gen++)
	{
//		fprintf(fgen, "%d   %6.4f+-%6.4f    %6.4f+-%6.4f    %6.4f+-%6.4f   %6.4f+-%6.4f   %6.4f+-%6.4f   %6.4f+-%6.4f   %6.4f+-%6.4f   %6.4f+-%6.4f   %6.4f+-%6.4f   %6.4f+-%6.4f\n",
//			gen, accmean(&gmean_s[gen]), se(&gmean_s[gen]), accmean(&gvar_s[gen]), se(&gvar_s[gen]), accmean(&SK2[gen]), se(&SK2[gen]), accmean(&QSrep[gen]), se(&QSrep[gen]), accmean(&HOrep[gen]), se(&HOrep[gen]), accmean(&HSrep[gen]), se(&HSrep[gen]), accmean(&HTrep[gen]), se(&HTrep[gen]), accmean(&FST[gen]), se(&FST[gen]), accmean(&NeH[gen]), se(&NeH[gen]), accmean(&NeFst[gen]), se(&NeFst[gen]));
        
		fprintf(fgen, "%d	%6.4f	%6.4f	%6.4f	%6.4f	%6.4f	%6.4f	%6.4f	%6.4f	%6.4f	%6.4f	%6.4f	%6.4f	%6.4f	%6.4f	%6.4f	%6.4f	%6.4f	%6.4f	%6.4f	%6.4f	%6.4f	%6.4f	%6.4f	%6.4f	%6.4f	%6.4f	%6.4f\n", gen, accmean(&gmean_s[gen]), accmean(&gmean_s_F[gen]), accmean(&ID_HOM[gen]), -accmean(&IDFyang2[gen]), accmean(&Wmax[gen]), (accmean(&Wmax[gen])-accmean(&gmean_s[gen]))/accmean(&Wmax[gen]), accmean(&gvar_s[gen]), accmean(&SK2[gen]), accmean(&AVE_VA[gen]), accmean(&AVE_VAeff[gen]), accmean(&AVE_VA_del[gen]), accmean(&AVE_VA_let[gen]), accmean(&AVE_VA_od[gen]),accmean(&AVE_VD[gen]), accmean(&AVE_VDeff[gen]), accmean(&AVE_VD_del[gen]), accmean(&AVE_VD_let[gen]), accmean(&AVE_VD_od[gen]), accmean(&AVE_ID[gen]), accmean(&AVE_B_FS[gen]), accmean(&AVE_ID_del[gen]), accmean(&AVE_ID_let[gen]), accmean(&AVE_ID_od[gen]), accmean(&FST[gen]), (4*NIND)/(2.0+accmean(&SK2[gen])), accmean(&NeH[gen]), accmean(&NeFst[gen]));  
	}
    
//	fprintf(fgen, "\n\n Number of crossovers = %6.4f\n", accmean(&numberofrecom));

	fprintf(fgen, "\ngen	VAdelr	VAodr     VDdelr    VDodr\n");
	for (gen=0; gen<=generations; gen++)
	{
		fprintf(fgen, "%d	%6.4f	%6.4f	%6.4f	%6.4f\n",
			gen, accmean(&AVE_VA_del_real[gen]), accmean(&AVE_VA_od_real[gen]), accmean(&AVE_VD_del_real[gen]), accmean(&AVE_VD_od_real[gen]));
	}

	fprintf(fgen,"\ngen	  n.del	   qdel	   n.let	   qlet	  n.od	  qod    fixsd	  lostsd\n");
	for (gen=0; gen<=generations; gen++)
	{
//		fprintf(fgen, "%d   %6.4f+-%6.4f    %6.4f+-%6.4f    %6.4f+-%6.4f   %6.4f+-%6.4f   %6.4f+-%6.4f   %6.4f+-%6.4f   %6.4f+-%6.4f   %6.4f+-%6.4f\n",
//			gen, accmean(&gmean_s[gen]), se(&gmean_s[gen]), accmean(&gvar_s[gen]), se(&gvar_s[gen]), accmean(&QSrep[gen]), se(&QSrep[gen]), accmean(&HSrep[gen]), se(&HSrep[gen]), accmean(&HTrep[gen]), se(&HTrep[gen]), accmean(&FST[gen]), se(&FST[gen]), accmean(&NeH[gen]), se(&NeH[gen]), accmean(&NeFst[gen]), se(&NeFst[gen]));
		fprintf(fgen, "%d    %6.4f   %6.4f   %6.4f   %6.4f   %6.4f   %6.4f   %6.4f   %6.4f\n",
		gen, accsum(&del[gen])/(replicates*lines), accmean(&qdel[gen]), accsum(&let[gen])/(replicates*lines), accmean(&qlet[gen]), accsum(&od[gen])/(replicates*lines), accmean(&qod[gen]), accsum(&fix_sd[gen])/(replicates*lines), accsum(&lost_sd[gen])/(replicates*lines));
	}

	fprintf(fgen,"\ngen	  sdelX	 sdelr	 hsdelX	hsdelr     BdelX     Bdelr   sod_r    Bo_X    Bo_r\n");
	for (gen=0; gen<=generations; gen++)
	{
        fprintf(fgen, "%d    %6.4f   %6.4f   %6.4f   %6.4f   %6.4f   %6.4f   %6.4f	%6.4f	%6.4f\n",
		gen, accmean(&sdel_X[gen]), accmean(&sdel_real[gen]), accmean(&hsdel_X[gen]), accmean(&hsdel_real[gen]), accmean(&AVE_ID_del_X[gen]), accmean(&AVE_ID_del_real[gen]), (accmean(&sod1_real[gen])+accmean(&sod2_real[gen]))/2.0, accmean(&AVE_ID_od_X[gen]), accmean(&AVE_ID_od_real[gen]));
	}

	fclose(fgen);
}

/* ***************************************************** */
