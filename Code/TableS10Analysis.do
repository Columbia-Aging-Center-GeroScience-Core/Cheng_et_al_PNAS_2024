


//Desktop
global OD "/Users/danielbelsky/Library/CloudStorage/OneDrive-cumc.columbia.edu"

//Origin and Destination files for datasets
	//DHWFS folder where Bertie uploaded data 
global origin "$OD/DHWFS/Data/Phenotype/HealthPhenotypes_ChengPNAS" 
	//Paper folder where Mengling 
global destination "$OD/Projects/MenglingCheng/DHWFS/Data/LumeyFiles"


foreach X in anthro CESD cognition Demographic diabetes ecg famhist health sf35a smoking{ 
	global X `X'
	import spss using "$origin/$X.sav", clear 
	save "$destination/$X.dta", replace 
	}
	
//Backbone Data 
import delimited using "$OD/DHWFS/Backbone/DHWFS_backbone.csv", delim(comma) varn(1) clear 
keep hofnum1 sex ageint agelab idfam subtype bproxdate almpdate  catthree3 anywk_new2 dwk_neg1vsall_new dwk_neg1vsall_new2 dwk0vsall_new2 dwk1vsall_new2 dwk2vsall_new2 dwk3vsall_new2 dwk4vsall_new2 prenatal_exp_days prenatal_exp_weeks prenatal_exp_timing_days prenatal_exp_timing_weeks postnatal_exp_days postnatal_exp_weeks preconceptional_exp_days preconceptional_exp_weeks gest_length

//Blood Pressure, Anthopometry, Hypertension 
merge 1:1 hofnum1 using "$destination/anthro.dta", nogen keepus(sbp dbp waist whr bmi)
recode sbp (0/139.999=0) (140/500=1), gen(hyp03)
replace hyp03 = 1 if dbp>=90 & dbp!=.
label var hyp03 "Hypertension 2003 Guidelines"
recode sbp (0/129.999=0) (130/500=1), gen(hyp17)
replace hyp17 = 1 if dbp>=80 & dbp!=.
label var hyp17 "Hypertension 2017 Guidelines"
recode bmi (0/29.9999=0) (30/500=1), gen(obese)
//Depression
merge 1:1 hofnum1 using "$destination/CESD.dta", nogen keepus(cesd CESDDX)
//Cognition 
merge 1:1 hofnum1 using "$destination/cognition.dta", nogen keepus(GENCOG)
//SES
merge 1:1 hofnum1 using "$destination/Demographic.dta", nogen keepus(pareduc educ5 nvrmarr employ)
//Diabetes
merge 1:1 hofnum1 using "$destination/diabetes.dta", nogen keepus(DIABETES)
//CVD
merge 1:1 hofnum1 using "$destination/ecg.dta", nogen keepus(CADQ FRAMRISK FRAMPT)
//FAMHIS
merge 1:1 hofnum1 using "$destination/famhist.dta", nogen keepus(FAMHIS)
//more CVD 
merge 1:1 hofnum1 using "$destination/health.dta", nogen keepus(mi stroke)
rename mi MI 
//Self-Rated Health
merge 1:1 hofnum1 using "$destination/sf35a.dta", nogen keepus(MCS PCS)
//Smoking 
merge 1:1 hofnum1 using "$destination/smoking.dta", nogen keepus(SMOKING PYTOTAL)

gen hofnum = hofnum1  
//Clocks
merge 1:1 hofnum using "$OD/Projects/MenglingCheng/DHWFS/Data/DNAmData/CPRDNAmClocks/CPRDNAmClocks_wZhang.dta",  keepus(age dnamage dnamagehannum dnamageskinblood dnamphenoage dnamgrimage pchorvath1 pchannum pcskinblood pcphenoage pcgrimage dunedinpace cd8t cd4t nk bcell mono gran propneuron plasmablast cd8pcd28ncd45ran cd8naive cd4naive idfam zhangclock)

drop if hofnum1==.

destring agelab, replace force 

foreach x in dnamage dnamagehannum dnamageskinblood dnamphenoage dnamgrimage zhangclock pchorvath1 pchannum pcskinblood pcphenoage pcgrimage dunedinpace {
	if `"`x'"' != "dunedinpace" { 
		reg `x' agelab
		capture drop r_`x' 
		predict r_`x', r
		capture drop z`x'
		egen z`x'=std(r_`x')	
		}
		capture drop z`x' 
		egen z`x'=std(`x')		
	}
recode bmi (0/24.999=0) (25/29.999=1) (30/100=2), gen(obcat)
capture drop cage 
gen cage = agelab-60
global cells "cd8t cd4t nk bcell mono gran propneuron plasmablast cd8pcd28ncd45ran cd8naive cd4naive"

global C1 "c.cage##c.cage sex"
global C2 "c.cage##c.cage sex $cells" 


//Prevalence of Chronic Diseases 
tabstat hyp03 hyp17 DIABETES MI stroke if dunedinpace!=., by(anywk_new2) save 
matrix CDsumstats = r(Stat1) \ r(Stat2) \ r(StatTotal)
matrix rownames CDsumstats = Control Famine AllDHWFS
matrix colnames CDsumstats = Hyp03 Hyp17 T2D MI Strk
matrix list CDsumstats 

putexcel set "/Users/danielbelsky/Library/CloudStorage/OneDrive-cumc.columbia.edu/Projects/MenglingCheng/DHWFS/MS/PNAS/R1/Final/R1Analysis.xlsx", modify
putexcel B2 = matrix(CDsumstats), names

//SRH Distribution
tabstat PCS if dunedinpace!=., by(anywk_new2) s(mean sd min max) save
matrix SRH = r(Stat1)' \ r(Stat2)' \ r(StatTotal)'
matrix rownames SRH = Control Famine AllDHWFS
matrix list SRH  

putexcel set "/Users/danielbelsky/Library/CloudStorage/OneDrive-cumc.columbia.edu/Projects/MenglingCheng/DHWFS/MS/PNAS/R1/Final/R1Analysis.xlsx", modify
putexcel B10 = matrix(SRH), names


//Compare famine effects on clocks with and without disease adjustment 
capture drop CDvars 
egen CDvars = rownonmiss(hyp03 DIABETES MI stroke)
xtset idfam 
foreach y in pcphenoage pcgrimage dunedinpace {
	//Famine Effect
	quietly xtgee z`y' anywk_new2 $C1 if dunedinpace!=. & CDvars==4, family(gaus) link(id) 
	matrix A = _b[anywk_new2], _b[anywk_new2]-1.96*_se[anywk_new2], _b[anywk_new2]+1.96*_se[anywk_new2], e(N), e(N_g)
	estimates store m1 
	quietly xtgee z`y' anywk_new2 $C1  hyp03 DIABETES MI stroke if dunedinpace!=., family(gaus) link(id) 
	matrix A = A \ _b[anywk_new2], _b[anywk_new2]-1.96*_se[anywk_new2], _b[anywk_new2]+1.96*_se[anywk_new2], e(N), e(N_g)	
	estimates store m2 
	esttab m1 m2, keep(anywk_new2) ci fixed compress nonum mtitle(Base Adj) title(`y') 
	matrix colnames A = b lb ub N Ng 
	matrix rownames A = base adj 
	matrix `y' = A 	
	}
	
matrix X = dunedinpace \ pcphenoage  \ pcgrimage 	
matrix rownames X = PACE adj Pheno adj Grim adj 
matrix list X 
putexcel set "/Users/danielbelsky/Library/CloudStorage/OneDrive-cumc.columbia.edu/Projects/MenglingCheng/DHWFS/MS/PNAS/R1/Final/R1Analysis.xlsx", modify 
putexcel B20 = matrix(X), names

//Self-rated Health Analysis
xtgee PCS anywk_new2 c.agelab##c.agelab sex if dunedinpace!=., family(gaus) link(id) 
matrix A = _b[anywk_new2], _b[anywk_new2]-1.96*_se[anywk_new2], _b[anywk_new2]+1.96*_se[anywk_new2], e(N), e(N_g)	
matrix rownames A = SRH
matrix colnames A = b lb ub N Ng 
matrix list A 
putexcel set "/Users/danielbelsky/Library/CloudStorage/OneDrive-cumc.columbia.edu/Projects/MenglingCheng/DHWFS/MS/PNAS/R1/Final/R1Analysis.xlsx", modify 
putexcel B40 = matrix(A)


//Compare famine effect on disease outcome w/ and w/o adj for DunedinPACE 
xtset idfam 
foreach y in hyp03 hyp17 DIABETES MI stroke {
	 xtgee `y' anywk_new2 $C1 if dunedinpace!=., family(binomial) link(logit) 
	matrix A = _b[anywk_new2], _b[anywk_new2]-1.96*_se[anywk_new2], _b[anywk_new2]+1.96*_se[anywk_new2], e(N), e(N_g)
	estimates store m1 
	quietly xtgee `y' anywk_new2 $C1 dunedinpace, family(binomial) link(logit)
	matrix A = A \ _b[anywk_new2], _b[anywk_new2]-1.96*_se[anywk_new2], _b[anywk_new2]+1.96*_se[anywk_new2], e(N), e(N_g)
	estimates store m2 
	esttab m1 m2, keep(anywk_new2) ci fixed compress nonum mtitle(Base Adj) title(`y') 
	matrix colnames A = b lb ub N Ng 
	matrix rownames A = base adj 
	matrix `y' = A 
	}
matrix Fx = hyp03 \ hyp17 \DIABETES \MI \stroke
matrix rownames Fx = Hyp03 Hyp17 T2D MI Strk  

