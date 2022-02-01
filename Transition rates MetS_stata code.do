************************
* 					   *
* Microsimulation MetS *
*				   	   *
************************

********************************
* This do-file consists of     *
* the codes to derive the      *
* transition rates for the     *
* microsimulation MetS paper   *
********************************

use "G:\OV19_0493\Paper Microsimulatie\Datasets\merged_MetS_1.dta", clear

* POPULATION CONSISTS OF ALL 18-65 YEAR OLDS

mi rename EDUCATION1 EDUCATION_1

net from "I:\Lifelines\Programs\STATA-packages\mimrgns"
net install mimrgns, replace

*************************************************
* 					STEP 1 						*
*************************************************
* Estimate the distribution of MetS in the starting population (e.g. at age 18) *
mi estimate, or level(95): logit METS_1new_wmed AGE i.SEX i.EDUCATION_1 if AGE<=27
* Calculate the prevalences for the different sexes and education categories *
mimrgns, at(AGE=18 SEX=1 EDUCATION_1=4) predict(pr) cmdmargins
mimrgns, at(AGE=18 SEX=2 EDUCATION_1=4) predict(pr) cmdmargins
mimrgns, at(AGE=18 SEX=1 EDUCATION_1=8) predict(pr) cmdmargins
mimrgns, at(AGE=18 SEX=2 EDUCATION_1=8) predict(pr) cmdmargins

*************************************************
* 					STEP 2 						*
*************************************************
* Estimate age-specific rates of entering MetS *

* keep only if no MetS at baseline
keep if METS_1new_wmed == 0

* distributions for MetS at second assessment and smoking at baseline
tab METS_4new_womedtia
tab SMOKING_1
tab BHLS_dich
tab ALCOHOL_1
tab DIET

* Estimate age-specific rates of entering MetS *
* by using agegroups
mi estimate, or level(95): logit METS_4new_womedtia i.AGEgroup i.SEX i.EDUCATION_1 i.SMOKING_1 i.BHLS_dich i.ALCOHOL_1 i.DIET INT_1A1_2A1_com

* Distribution smoking among low and high educated *
tab SMOKING_1 if EDUCATION_1 == 4
tab SMOKING_1 if EDUCATION_1 == 8

* Distribution health literacy among low and high educated *
tab BHLS_dich if EDUCATION_1 == 4
tab BHLS_dich if EDUCATION_1 == 8

* Distribution alcohol among low and high educated *
tab ALCOHOL_1 if EDUCATION_1 == 4
tab ALCOHOL_1 if EDUCATION_1 == 8

* Distribution diet among low and high educated *
tab DIET if EDUCATION_1 == 4
tab DIET if EDUCATION_1 == 8

* Generate smoking dummy's to make it easier to use for the margins *
mi passive: gen SMOKING_form = 0
mi passive: replace SMOKING_form = 1 if SMOKING_1 == 1
mi passive: replace SMOKING_form = . if SMOKING_1 == .
tab SMOKING_form SMOKING_1

mi passive: gen SMOKING_cur = 0
mi passive: replace SMOKING_cur = 1 if SMOKING_1 == 2
mi passive: replace SMOKING_cur = . if SMOKING_1 == .
tab SMOKING_cur SMOKING_1

* Generate health literacy dummy's to make it easier to use for the margins *
mi passive: gen HEALTHLIT_low = 0
mi passive: replace HEALTHLIT_low = 1 if BHLS_dich == 0
mi passive: replace HEALTHLIT_low = . if BHLS_dich == .
tab HEALTHLIT_low BHLS_dich

* Generate alcohol dummy's to make it easier to use for the margins *
mi passive: gen ALCOHOLUSE_no = 0
mi passive: replace ALCOHOLUSE_no = 1 if ALCOHOL_1 == 0
mi passive: replace ALCOHOLUSE_no = . if ALCOHOL_1 == .
tab ALCOHOLUSE_no ALCOHOL_1

mi passive: gen ALCOHOLUSE_high = 0
mi passive: replace ALCOHOLUSE_high = 1 if ALCOHOL_1 == 2
mi passive: replace ALCOHOLUSE_high = . if ALCOHOL_1 == .
tab ALCOHOLUSE_high ALCOHOL_1

* Generate diet dummy's to make it easier to use for the margins *
mi passive: gen DIET_low = 0
mi passive: replace DIET_low = 1 if DIET == 2
mi passive: replace DIET_low = . if DIET == .
tab DIET_low DIET

mi passive: gen DIET_moderate = 0
mi passive: replace DIET_moderate = 1 if DIET == 1
mi passive: replace DIET_moderate = . if DIET == .
tab DIET_moderate DIET

*********************
mi estimate, or level(95): logit METS_4new_womedtia i.AGEgroup i.SEX i.EDUCATION_1 SMOKING_form SMOKING_cur HEALTHLIT_low ALCOHOLUSE_no ALCOHOLUSE_high DIET_low DIET_moderate
* Calculate the incidences for the different sexes, education categories, and education-distributions for smoking health literacy alcohol use and diet*
mimrgns, at(AGEgroup=(1(1)10) SEX=1 EDUCATION_1=4 SMOKING_form=0.3792 SMOKING_cur=0.2384 HEALTHLIT_low=0.3194 ALCOHOLUSE_no=0.2235 ALCOHOLUSE_high=0.3611 DIET_low=0.1129 DIET_moderate=0.8040) predict(pr) cmdmargins
mimrgns, at(AGEgroup=(1(1)10) SEX=2 EDUCATION_1=4 SMOKING_form=0.3792 SMOKING_cur=0.2384 HEALTHLIT_low=0.3194 ALCOHOLUSE_no=0.2235 ALCOHOLUSE_high=0.3611 DIET_low=0.1129 DIET_moderate=0.8040) predict(pr) cmdmargins
mimrgns, at(AGEgroup=(1(1)10) SEX=1 EDUCATION_1=8 SMOKING_form=0.2350 SMOKING_cur=0.1143 HEALTHLIT_low=0.0835 ALCOHOLUSE_no=0.1142 ALCOHOLUSE_high=0.3988 DIET_low=0.0525 DIET_moderate=0.8152) predict(pr) cmdmargins
mimrgns, at(AGEgroup=(1(1)10) SEX=2 EDUCATION_1=8 SMOKING_form=0.2350 SMOKING_cur=0.1143 HEALTHLIT_low=0.0835 ALCOHOLUSE_no=0.1142 ALCOHOLUSE_high=0.3988 DIET_low=0.0525 DIET_moderate=0.8152) predict(pr) cmdmargins

* the before mentioned models are without follow-up time (INT_1A1_2A1_com) *
* calculate follow-up time *
gen time = INT_1A1_2A1_com/12
summarize INT_1A1_2A1_com time
tab time // time between baseline and second assessment = 3.96 year
* divide the estimate out of the model by 3,96 to get the MetS incidence per year

*************************************************************************
* step 4, estimate counterfactual age-specific rates of developing MetS *
*************************************************************************
* Based on the models estimated in step 2 and 3, calculate the age-specific rate of developing and remitting MetS for VMBO (EDUC=10) for *
* men and women if VMBO respondents would have the SAME SMOKING HABITS as university respondents *

* SMOKING *
mi estimate, or level(95): logit METS_4new_womedtia i.AGEgroup i.SEX i.EDUCATION_1 SMOKING_form SMOKING_cur HEALTHLIT_low ALCOHOLUSE_no ALCOHOLUSE_high DIET_low DIET_moderate
mimrgns, at(AGEgroup=(1(1)10) SEX=1 EDUCATION_1=4 SMOKING_form=0.2350 SMOKING_cur=0.1143 HEALTHLIT_low=0.3194 ALCOHOLUSE_no=0.2235 ALCOHOLUSE_high=0.3611 DIET_low=0.1129 DIET_moderate=0.8040) predict(pr) cmdmargins
mimrgns, at(AGEgroup=(1(1)10) SEX=2 EDUCATION_1=4 SMOKING_form=0.2350 SMOKING_cur=0.1143 HEALTHLIT_low=0.3194 ALCOHOLUSE_no=0.2235 ALCOHOLUSE_high=0.3611 DIET_low=0.1129 DIET_moderate=0.8040) predict(pr) cmdmargins

* HEALTH LITERACY *
mi estimate, or level(95): logit METS_4new_womedtia i.AGEgroup i.SEX i.EDUCATION_1 SMOKING_form SMOKING_cur HEALTHLIT_low ALCOHOLUSE_no ALCOHOLUSE_high DIET_low DIET_moderate
mimrgns, at(AGEgroup=(1(1)10) SEX=1 EDUCATION_1=4 SMOKING_form=0.3792 SMOKING_cur=0.2384 HEALTHLIT_low=0.0835 ALCOHOLUSE_no=0.2235 ALCOHOLUSE_high=0.3611 DIET_low=0.1129 DIET_moderate=0.8040) predict(pr) cmdmargins
mimrgns, at(AGEgroup=(1(1)10) SEX=2 EDUCATION_1=4 SMOKING_form=0.3792 SMOKING_cur=0.2384 HEALTHLIT_low=0.0835 ALCOHOLUSE_no=0.2235 ALCOHOLUSE_high=0.3611 DIET_low=0.1129 DIET_moderate=0.8040) predict(pr) cmdmargins

* ALCOHOL *
mi estimate, or level(95): logit METS_4new_womedtia i.AGEgroup i.SEX i.EDUCATION_1 SMOKING_form SMOKING_cur HEALTHLIT_low ALCOHOLUSE_no ALCOHOLUSE_high DIET_low DIET_moderate
mimrgns, at(AGEgroup=(1(1)10) SEX=1 EDUCATION_1=4 SMOKING_form=0.3792 SMOKING_cur=0.2384 HEALTHLIT_low=0.3194 ALCOHOLUSE_no=0.1142 ALCOHOLUSE_high=0.3988 DIET_low=0.1129 DIET_moderate=0.8040) predict(pr) cmdmargins
mimrgns, at(AGEgroup=(1(1)10) SEX=2 EDUCATION_1=4 SMOKING_form=0.3792 SMOKING_cur=0.2384 HEALTHLIT_low=0.3194 ALCOHOLUSE_no=0.1142 ALCOHOLUSE_high=0.3988 DIET_low=0.1129 DIET_moderate=0.8040) predict(pr) cmdmargins

* DIET *
mi estimate, or level(95): logit METS_4new_womedtia i.AGEgroup i.SEX i.EDUCATION_1 SMOKING_form SMOKING_cur HEALTHLIT_low ALCOHOLUSE_no ALCOHOLUSE_high DIET_low DIET_moderate
mimrgns, at(AGEgroup=(1(1)10) SEX=1 EDUCATION_1=4 SMOKING_form=0.3792 SMOKING_cur=0.2384 HEALTHLIT_low=0.3194 ALCOHOLUSE_no=0.2235 ALCOHOLUSE_high=0.3611 DIET_low=0.0525 DIET_moderate=0.8152) predict(pr) cmdmargins
mimrgns, at(AGEgroup=(1(1)10) SEX=2 EDUCATION_1=4 SMOKING_form=0.3792 SMOKING_cur=0.2384 HEALTHLIT_low=0.3194 ALCOHOLUSE_no=0.2235 ALCOHOLUSE_high=0.3611 DIET_low=0.0525 DIET_moderate=0.8152) predict(pr) cmdmargins

* SMOKING, HEALTH LITERACY, ALCOHOL AND DIET *
mi estimate, or level(95): logit METS_4new_womedtia i.AGEgroup i.SEX i.EDUCATION_1 SMOKING_form SMOKING_cur HEALTHLIT_low ALCOHOLUSE_no ALCOHOLUSE_high DIET_low DIET_moderate
mimrgns, at(AGEgroup=(1(1)10) SEX=1 EDUCATION_1=4 SMOKING_form=0.2350 SMOKING_cur=0.1143 HEALTHLIT_low=0.0835 ALCOHOLUSE_no=0.1142 ALCOHOLUSE_high=0.3988 DIET_low=0.0525 DIET_moderate=0.8152) predict(pr) cmdmargins
mimrgns, at(AGEgroup=(1(1)10) SEX=2 EDUCATION_1=4 SMOKING_form=0.2350 SMOKING_cur=0.1143 HEALTHLIT_low=0.0835 ALCOHOLUSE_no=0.1142 ALCOHOLUSE_high=0.3988 DIET_low=0.0525 DIET_moderate=0.8152) predict(pr) cmdmargins

*************************************************
* 					STEP 3 						*
*************************************************
* Estimate age-specific rates of remitting MetS *

use "G:\OV19_0493\Paper Microsimulatie\Datasets\merged_MetS_1.dta", clear
mi rename EDUCATION1 EDUCATION_1

keep if METS_1new_wmed == 1
tab METS_4new_womedtia

* Change reference category of MetS at T4
recode METS_4new_womedtia (1 = 0) (0 = 1)
label define METS_4new_change 0 "Yes" 1 "No"
label value METS_4new_womedtia METS_4new_change
tab METS_4new_womedtia,nolab

* Distribution smoking among low and high educated *
tab SMOKING_1 if EDUCATION_1 == 4
tab SMOKING_1 if EDUCATION_1 == 8

* Distribution health literacy among low and high educated *
tab BHLS_dich if EDUCATION_1 == 4
tab BHLS_dich if EDUCATION_1 == 8

* Distribution alcohol among low and high educated *
tab ALCOHOL_1 if EDUCATION_1 == 4
tab ALCOHOL_1 if EDUCATION_1 == 8

* Distribution diet among low and high educated *
tab DIET if EDUCATION_1 == 4
tab DIET if EDUCATION_1 == 8

* Generate smoking dummy's to make it easier to use for the margins *
mi passive: gen SMOKING_form = 0
mi passive: replace SMOKING_form = 1 if SMOKING_1 == 1
mi passive: replace SMOKING_form = . if SMOKING_1 == .
tab SMOKING_form SMOKING_1

mi passive: gen SMOKING_cur = 0
mi passive: replace SMOKING_cur = 1 if SMOKING_1 == 2
mi passive: replace SMOKING_cur = . if SMOKING_1 == .
tab SMOKING_cur SMOKING_1

* Generate health literacy dummy's to make it easier to use for the margins *
mi passive: gen HEALTHLIT_low = 0
mi passive: replace HEALTHLIT_low = 1 if BHLS_dich == 0
mi passive: replace HEALTHLIT_low = . if BHLS_dich == .
tab HEALTHLIT_low BHLS_dich

* Generate alcohol dummy's to make it easier to use for the margins *
mi passive: gen ALCOHOLUSE_no = 0
mi passive: replace ALCOHOLUSE_no = 1 if ALCOHOL_1 == 0
mi passive: replace ALCOHOLUSE_no = . if ALCOHOL_1 == .
tab ALCOHOLUSE_no ALCOHOL_1

mi passive: gen ALCOHOLUSE_high = 0
mi passive: replace ALCOHOLUSE_high = 1 if ALCOHOL_1 == 2
mi passive: replace ALCOHOLUSE_high = . if ALCOHOL_1 == .
tab ALCOHOLUSE_high ALCOHOL_1

* Generate diet dummy's to make it easier to use for the margins *
mi passive: gen DIET_low = 0
mi passive: replace DIET_low = 1 if DIET == 2
mi passive: replace DIET_low = . if DIET == .
tab DIET_low DIET

mi passive: gen DIET_moderate = 0
mi passive: replace DIET_moderate = 1 if DIET == 1
mi passive: replace DIET_moderate = . if DIET == .
tab DIET_moderate DIET

*************************************
mi estimate, or level(95): logit METS_4new_womedtia i.AGEgroup i.SEX i.EDUCATION_1 SMOKING_form SMOKING_cur HEALTHLIT_low ALCOHOLUSE_no ALCOHOLUSE_high DIET_low DIET_moderate
* Calculate the incidences for the different sexes, education categories, and education-distributions for smoking health literacy alcohol use and diet*
mimrgns, at(AGEgroup=(1(1)10) SEX=1 EDUCATION_1=4 SMOKING_form=0.4296 SMOKING_cur=0.2307 HEALTHLIT_low=0.3464 ALCOHOLUSE_no=0.2763 ALCOHOLUSE_high=0.3421 DIET_low=0.1121 DIET_moderate=0.7957) predict(pr) cmdmargins
mimrgns, at(AGEgroup=(1(1)10) SEX=2 EDUCATION_1=4 SMOKING_form=0.4296 SMOKING_cur=0.2307 HEALTHLIT_low=0.3464 ALCOHOLUSE_no=0.2763 ALCOHOLUSE_high=0.3421 DIET_low=0.1121 DIET_moderate=0.7957) predict(pr) cmdmargins
mimrgns, at(AGEgroup=(1(1)10) SEX=1 EDUCATION_1=8 SMOKING_form=0.4005 SMOKING_cur=0.1648 HEALTHLIT_low=0.1037 ALCOHOLUSE_no=0.1596 ALCOHOLUSE_high=0.4656 DIET_low=0.0660 DIET_moderate=0.8313) predict(pr) cmdmargins
mimrgns, at(AGEgroup=(1(1)10) SEX=2 EDUCATION_1=8 SMOKING_form=0.4005 SMOKING_cur=0.1648 HEALTHLIT_low=0.1037 ALCOHOLUSE_no=0.1596 ALCOHOLUSE_high=0.4656 DIET_low=0.0660 DIET_moderate=0.8313) predict(pr) cmdmargins

* the before mentioned models are without follow-up time (INT_1A1_2A1_com) *
* calculate follow-up time *
gen time = INT_1A1_2A1_com/12
summarize INT_1A1_2A1_com time
tab time // time between baseline and second assessment = 3.96 year
* divide the estimate out of the model by 3,96 to get the MetS incidence per year

*************************************************************************
* step 4, estimate counterfactual age-specific rates of developing MetS *
*************************************************************************
* Based on the models estimated in step 2 and 3, calculate the age-specific rate of developing and remitting MetS for VMBO (EDUC=10) for *
* men and women if VMBO respondents would have the SAME SMOKING HABITS as university respondents *

* SMOKING *
mi estimate, or level(95): logit METS_4new_womedtia i.AGEgroup i.SEX i.EDUCATION_1 SMOKING_form SMOKING_cur HEALTHLIT_low ALCOHOLUSE_no ALCOHOLUSE_high DIET_low DIET_moderate
mimrgns, at(AGEgroup=(1(1)10) SEX=1 EDUCATION_1=4 SMOKING_form=0.4005 SMOKING_cur=0.1648 HEALTHLIT_low=0.3464 ALCOHOLUSE_no=0.2763 ALCOHOLUSE_high=0.3421 DIET_low=0.1121 DIET_moderate=0.7957) predict(pr) cmdmargins
mimrgns, at(AGEgroup=(1(1)10) SEX=2 EDUCATION_1=4 SMOKING_form=0.4005 SMOKING_cur=0.1648 HEALTHLIT_low=0.3464 ALCOHOLUSE_no=0.2763 ALCOHOLUSE_high=0.3421 DIET_low=0.1121 DIET_moderate=0.7957) predict(pr) cmdmargins

* HEALTH LITERACY *
mi estimate, or level(95): logit METS_4new_womedtia i.AGEgroup i.SEX i.EDUCATION_1 SMOKING_form SMOKING_cur HEALTHLIT_low ALCOHOLUSE_no ALCOHOLUSE_high DIET_low DIET_moderate
mimrgns, at(AGEgroup=(1(1)10) SEX=1 EDUCATION_1=4 SMOKING_form=0.4296 SMOKING_cur=0.2307 HEALTHLIT_low=0.1037 ALCOHOLUSE_no=0.2763 ALCOHOLUSE_high=0.3421 DIET_low=0.1121 DIET_moderate=0.7957) predict(pr) cmdmargins
mimrgns, at(AGEgroup=(1(1)10) SEX=2 EDUCATION_1=4 SMOKING_form=0.4296 SMOKING_cur=0.2307 HEALTHLIT_low=0.1037 ALCOHOLUSE_no=0.2763 ALCOHOLUSE_high=0.3421 DIET_low=0.1121 DIET_moderate=0.7957) predict(pr) cmdmargins

* ALCOHOL *
mi estimate, or level(95): logit METS_4new_womedtia i.AGEgroup i.SEX i.EDUCATION_1 SMOKING_form SMOKING_cur HEALTHLIT_low ALCOHOLUSE_no ALCOHOLUSE_high DIET_low DIET_moderate
mimrgns, at(AGEgroup=(1(1)10) SEX=1 EDUCATION_1=4 SMOKING_form=0.4296 SMOKING_cur=0.2307 HEALTHLIT_low=0.3464 ALCOHOLUSE_no=0.1596 ALCOHOLUSE_high=0.4656 DIET_low=0.1121 DIET_moderate=0.7957) predict(pr) cmdmargins
mimrgns, at(AGEgroup=(1(1)10) SEX=2 EDUCATION_1=4 SMOKING_form=0.4296 SMOKING_cur=0.2307 HEALTHLIT_low=0.3464 ALCOHOLUSE_no=0.1596 ALCOHOLUSE_high=0.4656 DIET_low=0.1121 DIET_moderate=0.7957) predict(pr) cmdmargins

* DIET *
mi estimate, or level(95): logit METS_4new_womedtia i.AGEgroup i.SEX i.EDUCATION_1 SMOKING_form SMOKING_cur HEALTHLIT_low ALCOHOLUSE_no ALCOHOLUSE_high DIET_low DIET_moderate
mimrgns, at(AGEgroup=(1(1)10) SEX=1 EDUCATION_1=4 SMOKING_form=0.4296 SMOKING_cur=0.2307 HEALTHLIT_low=0.3464 ALCOHOLUSE_no=0.2763 ALCOHOLUSE_high=0.3421 DIET_low=0.0660 DIET_moderate=0.8313) predict(pr) cmdmargins
mimrgns, at(AGEgroup=(1(1)10) SEX=2 EDUCATION_1=4 SMOKING_form=0.4296 SMOKING_cur=0.2307 HEALTHLIT_low=0.3464 ALCOHOLUSE_no=0.2763 ALCOHOLUSE_high=0.3421 DIET_low=0.0660 DIET_moderate=0.8313) predict(pr) cmdmargins

* SMOKING, HEALTH LITERACY, ALCOHOL AND DIET *
mi estimate, or level(95): logit METS_4new_womedtia i.AGEgroup i.SEX i.EDUCATION_1 SMOKING_form SMOKING_cur HEALTHLIT_low ALCOHOLUSE_no ALCOHOLUSE_high DIET_low DIET_moderate
mimrgns, at(AGEgroup=(1(1)10) SEX=1 EDUCATION_1=4 SMOKING_form=0.4005 SMOKING_cur=0.1648 HEALTHLIT_low=0.1037 ALCOHOLUSE_no=0.1596 ALCOHOLUSE_high=0.4656 DIET_low=0.0660 DIET_moderate=0.8313) predict(pr) cmdmargins
mimrgns, at(AGEgroup=(1(1)10) SEX=2 EDUCATION_1=4 SMOKING_form=0.4005 SMOKING_cur=0.1648 HEALTHLIT_low=0.1037 ALCOHOLUSE_no=0.1596 ALCOHOLUSE_high=0.4656 DIET_low=0.0660 DIET_moderate=0.8313) predict(pr) cmdmargins
