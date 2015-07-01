#ifndef FITFUNC_HPP
#define FITFUNC_HPP

RooFitResult * fitTo(RooAbsPdf * pdf, RooDataSet * data, RooArgSet * cons = NULL, string opt = "" )
{
	//freopen ("tmp.txt","w",stdout);

	RooFitResult * res = NULL;
	if(cons!=NULL) 
		res = pdf->fitTo(*data,PrintLevel(-1),Save(),NumCPU(2),
			ExternalConstraints(*cons),Warnings(false), Verbose(kFALSE), Minos(true));
	else
		res = pdf->fitTo(*data,PrintLevel(-1),Save(),Warnings(false),NumCPU(2), Minos(true));
	
	//freopen ("/dev/tty", "a", stdout);
	return res;
}


void genPoints(unsigned np, double * a, double * f, double * ab)
{
	TRandom3 r(0);

	for(unsigned i=0; i<np; )
	{
		a[i] = r.Rndm()*1.5-0.75;
		f[i] = r.Rndm();
		ab[i] = r.Rndm() - 0.5;

		if( !((f[i]-1)*3./4. > a[i] || a[i] > -(f[i]-1)*3./4.) ) i++;
	}
}

void genPointsSphere(unsigned np, double * c, double *rad, double * a, double * f, double * ab)
{
	TRandom3 r(0);

	for(unsigned i=0; i<np; )
	{
		double myr = r.Rndm()*rad[0];
		double theta = r.Rndm()*2*TMath::Pi();
		double phi = r.Rndm()*TMath::Pi();
		a[i] = c[0] + myr*TMath::Cos(theta)*TMath::Sin(phi);
		f[i] = c[1] + myr*TMath::Sin(theta)*TMath::Sin(phi);
		ab[i] = c[2] + myr*TMath::Cos(phi);

		if( !((f[i]-1)*3./4. > a[i] || a[i] > -(f[i]-1)*3./4.) && TMath::Abs(ab[i]) <= 0.5 ) i++;
	}
}



RooFitResult * scan(RooAbsPdf * pdf, RooDataSet * data, RooAbsReal * nll, Str2VarMap p, int nfree, string opt, RooArgSet * cons, 
		unsigned np, double * a, double * f, double * ab, double *c, double *r, RooArgSet * nuisances)
{
	RooFitResult * res = NULL;
	double minLL = 1.e6, minLL2 = 1.e6;
	double mina2=0,minab2=0,minf2=0;

	bool afb_iscost  = p["afb"]->getAttribute("Constant");
	bool fL_iscost   = p["fL"]->getAttribute("Constant");
	bool afbB_iscost = p["afbB"]->getAttribute("Constant");
	((RooRealVar*)p["afb"])->setConstant();
	((RooRealVar*)p["fL"])->setConstant();
	((RooRealVar*)p["afbB"])->setConstant();

	for( unsigned i = 0; i < np; i++ )
	{
		if(!afb_iscost) ((RooRealVar*)p["afb"])->setVal(a[i]);
		if(!fL_iscost) ((RooRealVar*)p["fL"])->setVal(f[i]);
		if(!afbB_iscost) ((RooRealVar*)p["afbB"])->setVal(ab[i]);

		double tmp = nll->getVal();
		if(tmp < minLL)	{ minLL = tmp; c[0] = p["afb"]->getVal(); c[1] = p["fL"]->getVal(); c[2] = p["afbB"]->getVal(); }
		else if(tmp < minLL2) { minLL2 = tmp; mina2 = p["afb"]->getVal(); minf2 = p["fL"]->getVal(); minab2 = p["afbB"]->getVal(); }
	}

	if(!afb_iscost) ((RooRealVar*)p["afb"])->setVal(c[0]);
	if(!fL_iscost) ((RooRealVar*)p["fL"])->setVal(c[1]);
	if(!afbB_iscost) ((RooRealVar*)p["afbB"])->setVal(c[2]);
	if(nfree>3){ fixParam(pdf,new RooArgSet(),nuisances); res = fitTo(pdf,data,cons,opt); }
	((RooRealVar*)p["afb"])->setConstant(afb_iscost);
	((RooRealVar*)p["fL"])->setConstant(fL_iscost);
	((RooRealVar*)p["afbB"])->setConstant(afbB_iscost);

	c[0] = p["afb"]->getVal();
	c[1] = p["fL"]->getVal();
	c[2] = p["afbB"]->getVal();
	r[0] = TMath::Sqrt( TMath::Power(c[0]-mina2,2) + TMath::Power(c[1]-minf2,2) + TMath::Power(c[2]-minab2,2) );

	return res;
}



RooFitResult * findMin(RooAbsPdf * pdf, RooDataSet * data, RooAbsReal * nll, Str2VarMap p, vector< double > init, vector< double > end, int np, ISVALIDF_PTR isValid, int nfree = -1, string opt = "", RooArgSet * cons = NULL, RooArgSet * nuisances = NULL)
{
	RooFitResult * res = NULL;

	if(p.size()==2)
	{
		Str2VarMap::iterator it = p.begin(); 
		RooRealVar * p1 = (RooRealVar*)it->second;
		it++;
		RooRealVar * p2 = (RooRealVar*)it->second;

		double step = TMath::Abs(end[0]-init[0])/(double)np, step2 = TMath::Abs(end[1]-init[1])/(double)np;
		double minLogL = 1e6, min1 = 1e6, min2 = 1e6;

		bool afb_iscost  = p1->getAttribute("Constant");
		bool fL_iscost   = p2->getAttribute("Constant");
		p1->setConstant();
		p2->setConstant();

		if(!afb_iscost && !fL_iscost)
		{
			for(double a = init[0]; a <= end[0]; a+=step )
				for(double f = init[1]; f <= end[1]; f+=step2 )
				{
					if(!afb_iscost) p1->setVal(a);
					if(!fL_iscost) p2->setVal(f);

					if(isValid) if( !isValid(pdf) ) continue;

					double tmp = nll->getVal();
					if(tmp < minLogL) { minLogL = tmp; min1 = a; min2 = f; }	
				}

			p1->setVal(min1);
			p2->setVal(min2);
		}

		if(nfree>2 && opt.find("-nofit")==string::npos)
		{ 
			//fixParam(pdf,new RooArgSet(),nuisances); 
			res = fitTo(pdf,data,cons,opt);
		}
		p1->setConstant(afb_iscost);
		p2->setConstant(fL_iscost);
	}
	else if(p.size()==1)
	{
		double step = TMath::Abs(end[0]-init[0])/(double)np;
		double minLogL = 1e6, min1 = 1e6;

		Str2VarMap::iterator it = p.begin();
		RooRealVar * p1 = (RooRealVar*)it->second;
		bool iscost = p1->getAttribute("Constant");
		p1->setConstant();

		if(!iscost)
		{
			for(double a = init[0]; a <= end[0]; a+=step )
			{
				p1->setVal(a);

				if(isValid) if( !isValid(pdf) ) continue;
				double tmp = nll->getVal();
				if(tmp < minLogL) { minLogL = tmp; min1 = a; }
			}

			p1->setVal(min1);
		}

		if(nfree>1 && opt.find("-nofit")==string::npos)
		{
			//fixParam(pdf,new RooArgSet(),nuisances);
			res = fitTo(pdf,data,cons,opt);
		}
		p1->setConstant(iscost);
	}
	else cout << "ATTENTION: p size = " << p.size() << endl;

	return res;
}


RooFitResult * safeFit(RooAbsPdf * pdf, RooDataSet * data, Str2VarMap p, ISVALIDF_PTR isValid, string opt = "", int nfree = -1, RooArgSet * cons = NULL, RooAbsReal * nll = NULL)
{
	RooFitResult * res = NULL;

	RooRealVar cosThetaL("cosThetaL","cosThetaL",0.,-1.,1.);
	RooRealVar cosThetaB("cosThetaB","cosThetaB",0.,-1.,1.);
	
	RooArgSet obs(cosThetaL,cosThetaB);
	
	//if(opt.find("-scan")==string::npos) res = pdf->fitTo(*data,PrintLevel(-1),Save(),Extended(true)); 
	if(p.size()==1 && p.find("afb") != p.end())     p["fL"]  = GetParam(pdf,"fL");
	else if(p.size()==1 && p.find("fL") != p.end()) p["afb"] = GetParam(pdf,"afb");
	RooArgSet * nuisances = NULL;
	/*
	bool afb_iscost = false, fL_iscost = false, afbB_iscost = false;
	if (p.find("afb") != p.end())  { afb_iscost  = ((RooRealVar*)p["afb"])->getAttribute("Constant");  ((RooRealVar*)p["afb"])->setConstant();  }
	if (p.find("fL") != p.end())   { fL_iscost   = ((RooRealVar*)p["fL"])->getAttribute("Constant");   ((RooRealVar*)p["fL"])->setConstant();   }
	if (p.find("afbB") != p.end()) { afbB_iscost = ((RooRealVar*)p["afbB"])->getAttribute("Constant"); ((RooRealVar*)p["afbB"])->setConstant(); }
	RooArgSet * nuisances = copyFreePars(pdf,obs);
	if (p.find("afb") != p.end())  ((RooRealVar*)p["afb"])->setConstant(afb_iscost);
	if (p.find("afbB") != p.end()) ((RooRealVar*)p["afbB"])->setConstant(afbB_iscost);
	if (p.find("fL") != p.end())   ((RooRealVar*)p["fL"])->setConstant(fL_iscost);
	*/
	int np = 20;
	if((!res || res->covQual()!=3 || res->edm() > 0.1) && opt.find("-noscan")==string::npos)
	{
		if(!nll) nll = pdf->createNLL(*data);

		vector < double > mins, maxs, r;

		Str2VarMap::iterator iter; int pp = 0;
		for (iter = p.begin(); iter != p.end(); iter++) 
		{
			RooRealVar * curp = (RooRealVar *)iter->second;
			maxs.push_back(curp->getMax());
			mins.push_back(curp->getMin());
			r.push_back((maxs.back() - mins.back())/(double)np);
			pp++;
		}
		
		findMin(pdf,data,nll,p,mins,maxs,np,isValid,nfree,opt+"-nofit",cons,nuisances);
		
		double prec = 1e6;
		while (prec > 0.001)
		{
			double maxr = 0;
			maxs.clear(); mins.clear(); pp=0;
			for (iter = p.begin(); iter != p.end(); iter++) 
			{
				RooRealVar * curp = (RooRealVar *)iter->second;
				if((curp->getVal() + r[pp]) < curp->getMax()) maxs.push_back(curp->getVal() + r[pp]);
				else maxs.push_back(curp->getMax());
				if((curp->getVal() - r[pp]) > curp->getMin()) mins.push_back(curp->getVal() - r[pp]);
				else mins.push_back(curp->getMin());
				r[pp] = (maxs.back() - mins.back())/(double)np;
				if(r[pp] > maxr) maxr = r[pp];
				pp++;
			}
		
			prec = maxr;
			res = findMin(pdf,data,nll,p,mins,maxs,np,isValid,nfree,opt,cons,nuisances);
		}
		
		//if(!mynll) delete nll;
	}

	return res;
}


/*
RooFitResult * findMin(RooAbsPdf * pdf, RooDataSet * data, RooAbsReal * nll, Str2VarMap p, vector< double > init, vector< double > end, int np, ISVALIDF_PTR isValid, int nfree = -1, string opt = "", RooArgSet * cons = NULL)
{
	RooFitResult * res = NULL;

	if(p.size()==2)
	{
		double step = TMath::Abs(end[0]-init[0])/(double)np, step2 = TMath::Abs(end[1]-init[1])/(double)np;
		double minLogL = 1e6, min1 = 1e6, min2 = 1e6;

		bool afb_iscost  = p["afb"]->getAttribute("Constant");
		bool fL_iscost   = p["fL"]->getAttribute("Constant");
		((RooRealVar*)p["afb"])->setConstant();
		((RooRealVar*)p["fL"])->setConstant();

		if(!afb_iscost && !fL_iscost)
		{
			for(double a = init[0]; a <= end[0]; a+=step )
				for(double f = init[1]; f <= end[1]; f+=step2 )
				{
					((RooRealVar*)p["afb"])->setVal(a);
					((RooRealVar*)p["fL"])->setVal(f);

					if(isValid) if( !isValid(pdf) ) continue;

					double tmp = nll->getVal();
					if(tmp < minLogL) { minLogL = tmp; min1 = a; min2 = f; }	
				}

			((RooRealVar*)p["afb"])->setVal(min1);
			((RooRealVar*)p["fL"])->setVal(min2);
		}

		if(nfree>2 && opt.find("-nofit")==string::npos)
			res = fitTo(pdf,data,cons,opt);

		((RooRealVar*)p["afb"])->setConstant(afb_iscost);
		((RooRealVar*)p["fL"])->setConstant(fL_iscost);
	}
	else if(p.size()==1)
	{
		double step = TMath::Abs(end[0]-init[0])/(double)np;
		double minLogL = 1e6, min1 = 1e6;
		bool afbB_iscost = p["afbB"]->getAttribute("Constant");
		((RooRealVar*)p["afbB"])->setConstant();

		if(!afbB_iscost)
		{
			for(double a = init[0]; a <= end[0]; a+=step )
			{
				((RooRealVar*)p["afbB"])->setVal(a);

				if(isValid) if( !isValid(pdf) ) continue;
				double tmp = nll->getVal();
				if(tmp < minLogL) { minLogL = tmp; min1 = a; }
			}

			((RooRealVar*)p["afbB"])->setVal(min1);
		}

		if(nfree>1 && opt.find("-nofit")==string::npos) res = fitTo(pdf,data,cons,opt);
		((RooRealVar*)p["afbB"])->setConstant(afbB_iscost);
	}
	else cout << "ATTENTION: p size = " << p.size() << endl;

	return res;
}
*/


#endif
