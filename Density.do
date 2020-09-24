use "C:\Users\bluep\Dropbox\peng\Academia\Work with Brea\SNAD\SNAD data\Peng\SNAD-Analysis-T123-20200830-All imaging",clear
cd "C:\Users\bluep\Dropbox\peng\Academia\Work with Brea\SNAD\SNAD data\Peng\Density"
keep if time==1
recode ENSO (.=0),gen(enso)
egen cdrsum=rowtotal(cdrmem-cdrpers)
gen hippocampus=(Left_Hippocampus_t+Right_Hippocampus_t)/ICV_t
recode race (2/50=0),gen(white)

gen icv=ICV_t
gen amyg=Left_Amygdala_t + Right_Amygdala_t
gen hipo=Left_Hippocampus_t + Right_Hippocampus_t 
gen wmh=log(WM_hypointensities_t)
rename non_WM_hypointensities_t gmh
rename TotalGrayVol_t gmv
rename CerebralWhiteMatterVol_t wmv
gen medtemp=MeanMedTemp_Thick_old
gen frolobe=MeanFront_Thick_old
gen parlobe=MeanPar_Thick_old
gen temlobe=MeanTemp_Thick_old
gen occlobe=MeanOccip_Thick_old
reg MOCATOTS icv hipo amyg wmh frolobe parlobe temlobe occlobe ,vce(robust)
predict residual, residuals
egen reserve=std(residual) //no demographics

*compare density
label define density 0 "don't know" 1 "not very close" 2 "sort of close" 3 "very close"
label values density density
label var bdensity "Binary density of networks sort of close\very close"
label var b1density "Binary density of networks know each other"

preserve
keep if enso==0
tempfile h1 h2 h3
hist density, saving(`h1')  percent xlabel(0 1 2 3,valuelabel) 
hist bdensity, saving(`h2')  percent
hist b1density, saving(`h3')  percent
graph combine "`h1'" "`h2'" "`h3'"
graph export "Fig1a.jpg" , replace
restore

preserve
keep if enso==1
tempfile h1 h2 h3
hist density, saving(`h1')  percent xlabel(0 1 2 3,valuelabel) 
hist bdensity, saving(`h2')  percent
hist b1density, saving(`h3')  percent
graph combine "`h1'" "`h2'" "`h3'"
graph export "Fig1b.jpg" , replace
restore

bysort enso: pwcorr density bdensity b1density netsize MOCATOTS hippocampus reserve,sig

*regression
foreach x of varlist density bdensity b1density {
	egen `x'std=std(`x')
}
eststo clear
foreach x of varlist densitystd bdensitystd b1densitystd {
	eststo `x'1: reg MOCATOTS i.sex i.white age grade `x',vce(robust)
	eststo `x'2: reg hippocampus i.sex i.white age grade `x',vce(robust)
	eststo `x'3: reg reserve i.sex i.white age grade `x',vce(robust)
}
esttab densitystd1 bdensitystd1 b1densitystd1 ///
       densitystd2 bdensitystd2 b1densitystd2 /// 
	   densitystd3 bdensitystd3 b1densitystd3 /// 
  using "reg.csv", replace label nobaselevels ///
  b(%5.2f) se(%5.2f) star r2(%5.2f) aic(%5.2f) bic(%5.2f) nogap noconstant ///
  title("CRstd") noomitted 
