
#libraries
library(WDI)
library("ff")
library(stringr)
library(plyr)
library(tidyverse)
library(dplyr)
library(tidyr)
library(Hmisc)
library(dtplyr) 
library(data.table)
library(tidyfast)
library(glue)
library(lubridate)
library(readr)
library(janitor)
library(tools)

###set date pull
date_pull<-as_date("2023-01-10") #Date of data pull


#' Import demographic data
#' @param file.name Path of the demographics data file (CDISC format)
#' @param dtplyr.step Return the output as \code{dtplyr_step} to avoid unnecessary future calls to \code{as_tibble} or \code{as.data.table}
#' @import dplyr tibble 
#' @return Formatted demographic data as a tibble or \code{dtplyr_step}
#' @export import.demographic.data
#' 

dm <- fread(file="Partner_DM_2023-01-10.csv", skip = 0)  #Import the DM data file. Rename with appropriate dataset
colnames(dm) <- tolower(colnames(dm))

import.demographic.data <- function(file.name, dtplyr.step = FALSE){
  wdi_dat <- WDI(indicator = c("NY.GDP.PCAP.KD", "SP.DYN.LE00.IN", "SP.DYN.IMRT.IN"),
                 start = 2020, end = 2020, extra = TRUE)%>%
    filter(region != "Aggregates")%>%
    select("Alpha_3"=iso3c,income,region)
  country.lookup <- ISOcodes::ISO_3166_1 %>% as_tibble%>%
    mutate(Name=case_when(!is.na(Common_name)~Common_name,
                          Name=="Lao People's Democratic Republic"~"Lao PDR",
                          TRUE~Name))%>%dplyr::select(Alpha_3, Name)%>%left_join(wdi_dat)
  out <- dm %>%
    ###delete patients duplicates
    group_by(usubjid) %>%
    mutate(count=1)%>%
    mutate(n = sum(count)) %>%
    filter(n == 1) %>%
    ungroup()%>%
    rename(date_admit=rfstdtc)%>%
    as.data.frame()%>%
    mutate(age_d=case_when(ageu=="MONTHS"~12,
                           ageu=="YEARS" ~ 1,
                           ageu=="DAYS" ~ 365.25,
                           TRUE~ NA_real_))%>%
    mutate(age2=age/age_d)%>%
    select(-(age))%>%
    rename(age=age2)%>%
    mutate(age=replace(age,age<0,NA))%>%
    mutate(ethnic = iconv(ethnic, to ="ASCII//TRANSLIT") %>% tolower()) %>%
    mutate(ethnic = str_remove_all(ethnic, "\\s*\\([^)]*\\)")) %>%
    mutate(ethnic = str_replace_all(ethnic, " - ", "_")) %>%
    mutate(ethnic = str_replace_all(ethnic, "-", "_")) %>%
    mutate(ethnic = str_replace_all(ethnic, "/| / ", "_")) %>%
    mutate(ethnic = str_replace_all(ethnic, " ", "_")) %>%
    mutate(ethnic = str_replace_all(ethnic, ",", "_")) %>%
    mutate(ethnic = replace(ethnic, ethnic == "n_a" | ethnic == "na" | ethnic == "", NA))%>%
    #mutate(studyid=substr(usubjid,1, 7))%>%
    mutate(sex = case_when(sex == "M" ~ "Male",
                           sex == "F" ~ "Female",
                           TRUE ~ NA_character_))%>%
    # mutate(date_admit=replace(date_admit,date_admit >date_pull,NA))%>%
    select(usubjid, studyid, date_admit, age, sex, ethnic, country)
  site_id_country<-out%>%
    mutate(country = replace(country, country == "", NA)) %>%
    left_join(country.lookup, by = c("country" = "Alpha_3")) %>%
    select(-country) %>%
    rename(country = Name) %>%
    filter(!is.na(country))%>%
    distinct(usubjid, country, .keep_all =T)%>%
    select(usubjid, 'country_2'=country, income, region)
  out<-out%>%
    left_join(site_id_country)%>%
    mutate(country=country_2)%>%dplyr::select(-c(country_2))
  if(dtplyr.step){
    return(out)
  } else {
    return(out %>% as_tibble())
  }
}

imp_dm <- import.demographic.data(dm, dtplyr.step = FALSE)
save(imp_dm, file = "imp_dm.rda")


#' Import microb data
#' @param file.name Path of the microbio data file (CDISC format)
#' @param dtplyr.step Return the output as \code{dtplyr_step} to avoid unnecessary future calls to \code{as_tibble} or \code{as.data.table}
#' @import dplyr tibble 
#' @return Formatted demographic data as a tibble or \code{dtplyr_step}
#' @export import.microb.data

mb<-read.csv("Partner_MB_2023-01-10.csv") #Import the MB data file
colnames(mb) <- tolower(colnames(mb))

import.microb.data <- function(file.name, dtplyr.step = FALSE){
  
  mb<-mb
  
  detection<- mb%>%
    filter(mbtstdtl=="DETECTION")%>%
    filter(mbtestcd=="CRONAVIR"|mbtestcd=="SARSCOV2")%>%
    mutate(mbstresc = case_when(mbstresc == "NO" ~ "NEGATIVE",
                                mbstresc == "NEGATIVE" ~ "NEGATIVE",
                                mbstresc == "POSITIVE" ~ "POSITIVE",
                                TRUE ~ NA_character_)) %>%
    mutate(mbtestcd = paste0("cov_det_",mbtestcd)%>% tolower%>%str_replace_all(" ", "_")) %>%
    arrange(desc(mbstresc))%>%
    distinct(usubjid, mbtestcd, .keep_all =T)%>% 
    as.data.table() %>%
    pivot_wider(id_cols = usubjid, names_from = mbtestcd, values_from = mbstresc) %>%
    as.data.frame()
  
  identification<-mb%>%
    filter(mbtstdtl=="IDENTIFICATION")%>%
    distinct(usubjid, mbstresc, .keep_all =T)%>% 
    filter(mbstresc=="SEVERE ACUTE RESPIRATORY SYNDROME CORONAVIRUS 2"|
             mbstresc=="CORONAVIRIDAE")%>%
    mutate(mbstresc=replace(mbstresc,mbstresc=="SEVERE ACUTE RESPIRATORY SYNDROME CORONAVIRUS 2","SARSCOV2"))%>%
    mutate(mbstresc=replace(mbstresc,mbstresc=="SEVERE ACUTE RESPIRATORY SYNDROME-RELATED CORONAVIRUS","SARSCOV2"))%>%
    mutate(mbstresc=replace(mbstresc,mbstresc=="CORONAVIRIDAE","CRONAVIR"))%>%
    mutate(result="POSITIVE")%>%
    mutate(mbstresc = paste0("cov_id_",mbstresc)%>%
             tolower%>%
             str_replace_all(" ", "_")) %>%
    as.data.table() %>%
    pivot_wider(id_cols = usubjid, names_from = mbstresc, values_from = result) %>%
    as.data.frame()
  
  out<-full_join(detection,identification)%>%
    mutate(cov_det_id="NEGATIVE")%>%
    mutate(cov_det_id=case_when(cov_det_cronavir=="POSITIVE"|
                                  cov_det_sarscov2=="POSITIVE"|
                                  cov_id_cronavir=="POSITIVE"|
                                  cov_id_sarscov2=="POSITIVE"~
                                  "POSITIVE",
                                is.na(cov_det_cronavir)&
                                  is.na(cov_det_sarscov2)&
                                  is.na(cov_id_cronavir)&
                                  is.na(cov_id_sarscov2)~
                                  NA_character_,
                                TRUE~cov_det_id))
  
  if(dtplyr.step){
    return(out)
  } else {
    return(out %>% as_tibble())
  }
}

imp_mb <- import.microb.data(mb, dtplyr.step = FALSE)
save(imp_mb, file = "imp_mb.rda")

#' Process data on pregnancy (as comorbidity)
#' @param file.name Path of the dispositions data file (CDISC format)
#' @param dtplyr.step Return the output as \code{dtplyr_step} to avoid unnecessary future calls to \code{as_tibble} or \code{as.data.table}
#' @import dplyr tibble stringr
#' @return Formatted pregnancy data as a tibble or \code{dtplyr_step}
#' @export process.pregnancy.data

rp<- fread(file="Partner_RP_2023-01-10.csv", skip = 0) #Import the RP data file
colnames(rp) <- tolower(colnames(rp))

process.pregnancy.data <- function(file.name, dtplyr.step = FALSE){
  comorbid_pregnancy <- rp%>%
    filter(rptestcd=="PREGIND") %>%
    mutate(comorbid_pregnancy=rpstresc) %>%
    mutate(comorbid_pregnancy = case_when(comorbid_pregnancy == "Y" ~ TRUE,
                                          comorbid_pregnancy == "N" ~ FALSE,
                                          TRUE ~ NA)) %>%
    select(usubjid,comorbid_pregnancy)
  if(dtplyr.step){
    return(comorbid_pregnancy %>% lazy_dt(immutable = FALSE))
  } else {
    return(comorbid_pregnancy %>% as_tibble())
  }
}

imp_rp <- process.pregnancy.data(rp, dtplyr.step = FALSE)
save(imp_rp, file = "imp_rp.rda")

#' Process data on vital sign
#' @param file.name Path of the dispositions data file (CDISC format)
#' @param dtplyr.step Return the output as \code{dtplyr_step} to avoid unnecessary future calls to \code{as_tibble} or \code{as.data.table}
#' @import dplyr tibble dtplyr tidyfast
#' @importFrom data.table as.data.table
#' @importFrom glue glue
#' @return Formatted vital sign (wide format) as a tibble or \code{dtplyr_step}
#' @export process.vital.sign.data

vs<-read.csv("Partner_VS_2023-01-10.csv") #import the VS data file
colnames(vs) <- tolower(colnames(vs))

process.vital.sign.data <- function(file.name, dtplyr.step = FALSE){
  vital_sign <- vs %>%
    select(usubjid, vstestcd, vscat,vsstresn,vsstresu, vso2src, vsdy) %>%
    filter(vscat=="SIGNS AND SYMPTOMS AT HOSPITAL ADMISSION" | vscat=="SIGNS AND SYMPTOMS AT ADMISSION")%>%
    
    mutate(vsstresn=as.numeric(vsstresn))%>%
    mutate(vsstresn=case_when(vstestcd=="OXYSAT"& vsstresn< 1~ NA_real_,
                              vstestcd=="OXYSAT"& vsstresn> 100~ NA_real_,
                              
                              vstestcd=="BMI"& vsstresn< 0~ NA_real_,
                              vstestcd=="BMI"& vsstresn> 100~ NA_real_,
                              
                              vstestcd=="DIABP"& vsstresn< 0~ NA_real_,
                              vstestcd=="DIABP"& vsstresn> 300~ NA_real_,
                              
                              vstestcd=="HEIGHT"& vsstresn< 0~ NA_real_,
                              vstestcd=="HEIGHT"& vsstresn> 250~ NA_real_,
                              
                              vstestcd=="HR"& vsstresn< 0~ NA_real_,
                              vstestcd=="HR"& vsstresn> 250~ NA_real_,
                              
                              vstestcd=="MAP"& vsstresn< 0~ NA_real_,
                              vstestcd=="MAP"& vsstresn> 250~ NA_real_,
                              
                              vstestcd=="MUARMCIR"& vsstresn< 0~ NA_real_,
                              vstestcd=="MUARMCIR"& vsstresn> 100~ NA_real_,
                              
                              vstestcd=="PULSE"& vsstresn< 0~ NA_real_,
                              vstestcd=="PULSE"& vsstresn> 250~ NA_real_,
                              
                              vstestcd=="RESP"& vsstresn< 0~ NA_real_,
                              vstestcd=="RESP"& vsstresn> 60~ NA_real_,
                              
                              vstestcd=="SYSBP"& vsstresn< 0~ NA_real_,
                              vstestcd=="SYSBP"& vsstresn> 250~ NA_real_,
                              
                              vstestcd=="TEMP"& vsstresn< 30~ NA_real_,
                              vstestcd=="TEMP"& vsstresn> 44~ NA_real_,
                              
                              vstestcd=="WEIGHT"& vsstresn< 0~ NA_real_,
                              vstestcd=="WEIGHT"& vsstresn> 300~ NA_real_,
                              
                              TRUE~vsstresn))%>%
    filter(!is.na(vsstresn))%>%
    arrange(desc(vsdy))%>%
    distinct(usubjid,vstestcd, .keep_all =T)%>%
    mutate(vso2src=case_when(vso2src==""&vstestcd=="OXYSAT"~'UNKNOWN',
                             TRUE~vso2src))%>%
    mutate(vso2src= str_replace_all(vso2src, " ", "_"))%>%
    mutate(vstestcd=case_when(vstestcd=="OXYSAT"~paste0(vstestcd,"_",vso2src),
                              TRUE~vstestcd))%>%
    mutate(vstestcd = paste0("vs_",vstestcd)) %>%
    mutate(vstestcd = iconv(vstestcd, to ="ASCII//TRANSLIT") %>% tolower()) %>%
    as.data.table() %>%
    pivot_wider(id_cols = usubjid, names_from = vstestcd,  values_from = vsstresn)%>%
    as.data.frame() %>%
    mutate(vs_oxysat=case_when(!is.na(vs_oxysat_oxygen_therapy)~vs_oxysat_oxygen_therapy,
                               !is.na(vs_oxysat_room_air)~vs_oxysat_room_air,
                               TRUE~vs_oxysat_unknown))
  
  if(dtplyr.step){
    return(vital_sign)
  } else {
    return(vital_sign %>% as_tibble())
  }
  
} 

imp_vs<- process.vital.sign.data(vs, dtplyr.step = FALSE)
save(imp_vs, file = "imp_vs.rda")

#' Process data on laboratory
#' @param file.name Path of the dispositions data file (CDISC format)
#' @param dtplyr.step Return the output as \code{dtplyr_step} to avoid unnecessary future calls to \code{as_tibble} or \code{as.data.table}
#' @import dplyr tibble dtplyr tidyfast
#' @importFrom data.table as.data.table
#' @importFrom glue glue
#' @return Formatted laboratory (wide format) as a tibble or \code{dtplyr_step}
#' @export process.laboratory.data

lb<-read.csv("Partner_LB_2023-01-10.csv") #import the LB data file
colnames(lb) <- tolower(colnames(lb))

process.laboratory.data <- function(file.name, dtplyr.step = FALSE){
  laboratory <- lb%>%
    dplyr::select(usubjid, lbdy, lbtestcd, lbcat,lborres) %>%
    mutate(lborres=replace(lborres,lborres=="",NA))%>%
    mutate(studyid=substr(usubjid,1, 7))%>%
    mutate(lbcat=case_when(lbdy==1 & (studyid=="CVCCPUK"| 
                                        studyid=="CVMEWUS"| 
                                        studyid=="CORE"|
                                        studyid=="CVTDWXD"|
                                        studyid=="CVTTYLU"|
                                        studyid=="CVZXZMV"|
                                        studyid=="CVKBQEI"|
                                        studyid=="CVVECMO"|
                                        studyid=="CVBCFGF") ~"LABORATORY RESULTS ON ADMISSION",
                           #lbdy==1 & studyid=="CVMEWUS"~"LABORATORY RESULTS ON ADMISSION",
                           TRUE~as.character(lbcat)))%>%
    filter(lbcat=="LABORATORY RESULTS ON ADMISSION" | lbcat=="LABORATORY RESULTS ON ADMISSION")%>%
    filter(lbtestcd=="ALT"|
             lbtestcd=="APTT"|
             lbtestcd=="CRP"|
             lbtestcd=="LYM"|
             lbtestcd=="NEUT"|
             lbtestcd=="PT"|
             lbtestcd=="WBC"|
             lbtestcd=="BILI"|
             lbtestcd=="AST"|
             lbtestcd=="UREAN")%>%
    mutate(lborres=as.numeric(lborres))%>%
    filter(!is.na(lborres))%>%
    distinct(usubjid,lbtestcd, .keep_all =T)%>%
    mutate(lborres=case_when(lbtestcd=="NEUT" & lborres>100 ~ lborres/1000,
                             
                             lbtestcd=="LYM" & lborres>100 ~ lborres/1000,
                             
                             lbtestcd=="WBC" & lborres>100 ~ lborres/1000, 
                             
                             lbtestcd=="ALT" & lborres>2000 ~ NA_real_,
                             lbtestcd=="ALT" & lborres<0 ~ NA_real_,
                             
                             lbtestcd=="AST" & lborres>2000 ~ NA_real_,
                             lbtestcd=="AST" & lborres<0 ~ NA_real_,
                             
                             lbtestcd=="BILI" & lborres>2000 ~ NA_real_,
                             lbtestcd=="BILI" & lborres<0 ~ NA_real_,
                             
                             lbtestcd=="CRP" & lborres>500 ~ NA_real_,
                             lbtestcd=="CRP" & lborres<0 ~ NA_real_,
                             
                             lbtestcd=="PT" & lborres>105 ~ NA_real_,
                             lbtestcd=="PT" & lborres<0 ~ NA_real_,
                             
                             lbtestcd=="UREAN" & lborres>100 ~ NA_real_,
                             lbtestcd=="UREAN" & lborres<0 ~ NA_real_,
                             
                             lbtestcd=="APTT" & lborres>2000 ~ NA_real_,
                             lbtestcd=="APTT" & lborres<0 ~ NA_real_,
                             
                             TRUE ~ lborres ))%>%
    mutate(lborres=case_when(lbtestcd=="NEUT" & lborres>100 ~ NA_real_,
                             lbtestcd=="NEUT" & lborres<0 ~ NA_real_,
                             
                             lbtestcd=="LYM" & lborres>100 ~ NA_real_,
                             lbtestcd=="LYM" & lborres<0 ~ NA_real_,
                             
                             lbtestcd=="WBC" & lborres>100 ~ NA_real_,
                             lbtestcd=="WBC" & lborres<0 ~ NA_real_,
                             
                             TRUE ~ lborres ))%>%
    
    mutate(lbtestcd  = paste0("lab_",lbtestcd )) %>%
    #mutate(lbtestcd = glue("lab_{lbtestcd}", lbtestcd = lbtestcd)) %>%
    mutate(lbtestcd = iconv(lbtestcd, to ="ASCII//TRANSLIT") %>% tolower()) %>%
    as.data.table() %>%
    pivot_wider(id_cols = usubjid, names_from = lbtestcd,  values_from = lborres)
  
  
  if(dtplyr.step){
    return(laboratory)
  } else {
    return(laboratory%>% as_tibble())
  }
  
}  

imp_lb <- process.laboratory.data(lb, dtplyr.step = FALSE)
save(imp_lb, file = "imp_lb.rda")

#' Process data on outcomes
#' @param file.name Path of the dispositions data file (CDISC format)
#' @param dtplyr.step Return the output as \code{dtplyr_step} to avoid unnecessary future calls to \code{as_tibble} or \code{as.data.table}
#' @import dplyr tibble stringr
#' @return Formatted outcome data (long format) as a tibble or \code{dtplyr_step}
#' @export process.outcome.data

ds<-read.csv("Partner_DS_2023-01-10.csv") #Import the DS data file
colnames(ds) <- tolower(colnames(ds))

process.outcome.data <- function(file.name, dtplyr.step = FALSE){
  outcome <- ds%>%
    select(usubjid, dsterm, dscdstdy, dsdecod, dsdy, dsstdy, dsmodify) %>%
    mutate(outcome=tolower(dsterm))%>%
    mutate(outcome=case_when(outcome=="palliative"~"transferred",
                             outcome=="transferred to another unit"~"ongoing care",
                             outcome=="Ongoing health care needs NOT related to COVID episode"~"discharge",
                             outcome==""~NA_character_,
                             TRUE~outcome))%>%
    mutate(outcome=case_when(outcome%like%"hospitalis"~"ongoing care",
                             outcome%like%"hospitaliz"~"ongoing care",
                             outcome%like%"ongoing"~"ongoing care",
                             outcome=="in hospital"~"ongoing care",
                             
                             outcome%like%"death"~"death",
                             outcome=="died"~"death",
                             outcome=="deceased"~"death",
                             outcome=="died (non-covid)"~"death",
                             
                             #outcome=="Death In Hospital"~"Death",
                             outcome=="alive"~"discharge",
                             outcome%like%"discharge"~"discharge",
                             outcome%like%"transfer"~"transferred",
                             outcome=="long term care facility"~"transferred",
                             outcome=="quarantine center"~"transferred",
                             outcome=="missing in database"~"unknown outcome",
                             outcome=="unknown"~"unknown outcome",
                             outcome=="not recorded"~"unknown outcome",
                             TRUE ~ outcome))%>%
    group_by(usubjid) %>% 
    mutate(count=1)%>% 
    mutate(n = sum(count)) %>%
    filter(n == 1)%>%
    select(-c(dsterm,dsmodify,n,count))
  
  
  if(dtplyr.step){
    return(outcome)
  } else {
    return(outcome %>% as_tibble())
  }
  
}

imp_ds <-process.outcome.data(ds, dtplyr.step = FALSE)
save(imp_ds, file = "imp_ds.rda")

#' Process data on ICU admission
#' @param file.name Path of the healthcare encounters data file (CDISC format)
#' @param dtplyr.step Return the output as \code{dtplyr_step} to avoid unnecessary future calls to \code{as_tibble} or \code{as.data.table}
#' @import dplyr tibble tidyfast dtplyr
#' @importFrom data.table as.data.table
#' @return Formatted symptom data as a tibble or \code{dtplyr_step}
#' @export process.ICU.data

ho<-read.csv("Partner_HO_2023-01-10.csv") #Import HO data file
colnames(ho) <- tolower (colnames(ho))

process.ICU.data <- function(file.name, dtplyr.step = FALSE){
  icu <- ho%>%
    mutate(hooccur = case_when(hooccur == "Y" ~ TRUE,
                               hooccur == "N" ~ FALSE,
                               TRUE ~ NA)) %>%
    filter(!is.na (hooccur))%>%
    select(usubjid, hodecod, hooccur, hody, hostdy,hoendy, hodur, hocdstdy)
  
  icu <-icu%>%
    filter(hodecod=="INTENSIVE CARE UNIT")%>%
    arrange(desc(hody))%>%
    distinct(usubjid, .keep_all =T)%>%
    rename(ever_icu=hooccur)%>%
    select(-c(hodecod))
  
  if(dtplyr.step){
    return(icu)
  } else {
    return(icu %>% as_tibble())
  }
}

imp_icu<- process.ICU.data(ho, dtplyr.step = FALSE)
save(imp_icu, file = "imp_icu.rda")

#' Import data on symptoms and comorbidities
#' @param file.name Path of the symptoms data file (CDISC format)
#' @param dtplyr.step Return the output as \code{dtplyr_step} to avoid unnecessary future calls to \code{as_tibble} or \code{as.data.table}
#' @import dplyr tibble stringr
#' @return Formatted comorbidity and symptom data as a tibble or \code{dtplyr_step}
#' @export import.symptom.and.comorbidity.data

sa <- fread(file="Partner_SA_v2_20230110.csv", skip = 0) #Import SA data file
sa<-as_tibble(sa)
colnames(sa) <- tolower (colnames(sa))

import.symptom.and.comorbidity.data <- function(file.name, dtplyr.step = TRUE){
  
  out<-sa%>%
    select(usubjid, studyid, saterm, sacat, samodify, sapresp, saoccur, sady, sastdy, sadur) %>%
    mutate(sacat=case_when(
      saterm=="CLINICALLY-DIAGNOSED COVID-19"~"CLINICALLY-DIAGNOSED COVID-19",
      TRUE~sacat))%>%
    filter(sacat=="MEDICAL HISTORY"|
             sacat=="SIGNS AND SYMPTOMS AT HOSPITAL ADMISSION"|
             sacat=="CLINICALLY-DIAGNOSED COVID-19"|
             sacat=="SIGNS AND SYMPTOMS AT INITIAL ACUTE COVID-19 ILLNESS")%>% #follow-up patients with acute data
    mutate(saoccur = case_when(saoccur == "Y" ~ TRUE,
                               saoccur == "N" ~ FALSE,
                               TRUE ~ NA)) %>%
    filter(!is.na(saoccur)) %>%
    mutate(saterm=toupper(saterm))%>%
    mutate(saterm=case_when(samodify!=""|is.na(samodify)~samodify,
                            TRUE ~ saterm))%>%
    
    mutate(saterm=case_when(saterm%like%'CARDIAC ARRHYTHMIA'~'CHRONIC CARDIAC DISEASE',
                            saterm%like%'CARDIAC DISEASE'~'CHRONIC CARDIAC DISEASE',
                            saterm%like%'CHORNIC CARDIAC DISEASE'~'CHRONIC CARDIAC DISEASE',
                            saterm%like%'CHRONIC HEART DISEASE'~'CHRONIC CARDIAC DISEASE',
                            saterm%like%'CONGENITAL CA'~'CHRONIC CARDIAC DISEASE',
                            saterm%like%'CONGENTIAL CARDIOPATHY'~'CHRONIC CARDIAC DISEASE',
                            saterm=='CORONARY DISEASE'~'CHRONIC CARDIAC DISEASE',
                            saterm=='HEART FAILURE'~'CHRONIC CARDIAC DISEASE',
                            saterm=='OROVALVA DISEASE'~'CHRONIC CARDIAC DISEASE',
                            saterm=='RHEUMATIC HEART DISEASE'~'CHRONIC CARDIAC DISEASE',
                            saterm=='VALVULAR HEART DISEASE'~'CHRONIC CARDIAC DISEASE',
                            saterm=='CONGESTIVE HEART FAILURE'~'CHRONIC CARDIAC DISEASE',
                            saterm=='CORONARY ARTERY DISEASE'~'CHRONIC CARDIAC DISEASE',
                            TRUE~saterm))%>%
    
    mutate(saterm=case_when(saterm=='CHRONIC DIALYSIS'~'CHRONIC KIDNEY DISEASE',
                            saterm%like%'DEPRESSION'~'PSYCHIATRIC CONDITION',
                            saterm%like%'PSYCHOSIS'~'PSYCHIATRIC CONDITION',
                            saterm%like%'DYSLIPIDEMIA'~'CHRONIC METABOLIC DISORDER',
                            saterm%like%'HYPOTHYROIDISM'~'CHRONIC ENDOCRINE DISORDER NON DIABITES',
                            saterm%like%'HEPATITIS'~'LIVER DISEASE',	
                            saterm%like%'MARASUMAS'~'MALNUTRITION',	
                            saterm=='SAM UNDEFINED'~'MALNUTRITION',	
                            saterm=='MIXED MARASMIC-KWASH'~'MALNUTRITION',	
                            saterm=='OSA/ HOME CPAP/BI-PAP USE'~'OBESITY',
                            saterm=='PAOD'~'OTHER COMORBIDITIES',
                            saterm=='PEPTIC ULCER DISEASE EXCLUDING BLEEDING'~'OTHER COMORBIDITIES',
                            saterm=='PARALYSIS'~'CHRONIC NEUROLOGICAL DISORDER',
                            saterm=='STROKE OR OTHER NEUROLOGICAL DISORDERS'~'CHRONIC NEUROLOGICAL DISORDER',
                            saterm=='PULMONARY CIRCULATION DISORDER'~'CHRONIC CARDIAC DISEASE',
                            saterm%like%'ARRHYTHMIA'~'CHRONIC CARDIAC DISEASE',
                            TRUE ~ saterm))%>%
    
    mutate(saterm=case_when(saterm=='SUBSTANCE USE DISORDER'~'SUBSTANCE MISUSE',	
                            saterm=='VENOUS THROMBOEMBOLISM- DVT/PE'~'THROMBOLIC DISORDERS',
                            
                            saterm=='CHILLS/RIGORS'~'RIGOR OR SWEATING',
                            saterm=='NIGHT SWEAT'~'RIGOR OR SWEATING',
                            
                            saterm=='CONGESTION/RHINORRHEA'~'RUNNY NOSE',
                            saterm=='CONJUNCTIVAL CONGESTION'~'UPPER RESPIRATORY TRACT SYMPTOMS',
                            saterm=='SNEEZING'~'UPPER RESPIRATORY TRACT SYMPTOMS',
                            
                            saterm=='DELIRIUM / ENCEPHALOPATHY'~'ALTERED CONSCIOUSNESS CONFUSION',
                            saterm=='DIZZINESS/LIGHTHEADEDNESS'~'OTHER SIGNS AND SYMPTOMS',	
                            saterm=='GASTROGASTROINTESTINAL HEMORRHAGE'~'OTHER SIGNS AND SYMPTOMS',	
                            TRUE~saterm))%>%
    
    mutate(saterm=case_when(saterm%like%'TUBERCULOSIS'~'TUBERCULOSIS',
                            saterm%like%'MALIGNANCY'~'MALIGNANT NEOPLASM',
                            saterm%like%'SPECIFIC CANCERS'~'MALIGNANT NEOPLASM',
                            saterm%like%'SOLID TUMOR'~'MALIGNANT NEOPLASM',
                            saterm%like%'METASTATIC CANCER'~'MALIGNANT NEOPLASM',
                            
                            saterm=='SORE THROAT/THROAT PAIN'~'SORE THROAT',
                            saterm=='COAGULOPATHY'~'CHRONIC HEMATOLOGIC DISEASE',
                            saterm=='DYSLIPIDEMIA/HYPERLIPIDEMIA'~'CHRONIC HEMATOLOGIC DISEASE',
                            saterm=='IRON DEFICIENCY ANEMIA'~'CHRONIC HEMATOLOGIC DISEASE',
                            saterm=='BLOOD LOSS ANEMIA'~'CHRONIC HEMATOLOGIC DISEASE',
                            
                            saterm=='CHRONIC HEMATOLOGICAL DISEASE'~'CHRONIC HEMATOLOGIC DISEASE',
                            saterm=='CHRONIC LIVER DISEASE'~'LIVER DISEASE',
                            saterm%like%'ACUTE LIVER'~'LIVER DISEASE',
                            saterm%like%'CHRONIC RENAL FAILURE'~'CHRONIC KIDNEY DISEASE',
                            saterm%like%'CHRONIC LUNG DISEASE'~'CHRONIC PULMONARY DISEASE',
                            saterm%like%'CHROMIC PULMONARY DISEASE'~'CHRONIC PULMONARY DISEASE',
                            TRUE~saterm))%>%
    
    mutate(saterm=case_when(saterm%like%'RHEMATOLOGICAL DISORDER'~'rheumatologic disorder',
                            saterm%like%'CHRONIC NEUROLOGICAL'~'CHRONIC NEUROLOGICAL DISORDER',
                            saterm%like%'CURRENT SMOK'~'SMOKING',
                            saterm%like%'DIABETES'~'DIABETES',
                            saterm=='HISTORY OF PERIPHERAL OR CARDIAC REVASCULARIZATION'~'HISTORY OF PERIPHERAL OR CARDIAC REVASCULARIZATION',
                            saterm=='HISTORY OF SMOKING'~'SMOKING',
                            saterm%like%'SMOKING'~'SMOKING',
                            saterm%like%'HIV'~'AIDS/HIV',
                            saterm%like%'LIVER DISEASE'~'LIVER DISEASE',
                            saterm%like%'OTHER RELEVANT RISK'~'OTHER COMORBIDITIES',
                            saterm=='OTHER RISK FACTOR'~'OTHER COMORBIDITIES',
                            saterm%like%'RHEUMATOLOGICAL DISORD'~'RHEUMATOLOGIC DISORDER',
                            saterm=='SMOKER'~'SMOKING',
                            saterm=='SMOKER - CURRENT'~'SMOKING',
                            saterm=='SMOKER - FORMER'~'SMOKING - FORMER',
                            saterm=='FEEDING INTOLERANCE (PAEDIATRICS)'~'ANOREXIA',
                            saterm=='REFUSING TO EAT OR DRINK/HISTORY OF POOR ORAL INTAKE'~'ANOREXIA',
                            saterm%like%'ANOREXIA'~'ANOREXIA',
                            saterm=='ANOREXIA - LOSS OF APPETITE'~'ANOREXIA',
                            saterm=='CHEST PAIN/TIGHTNESS'~'CHEST PAIN',
                            saterm=='SWOLLEN NECK GLANDS/LYMPHADENOPATHY'~'LYMPHADENOPATHY',
                            TRUE~saterm))%>%
    
    mutate(saterm=case_when(saterm%like%'COUGH'~'COUGH',
                            saterm%like%'COUTH'~'COUGH',
                            saterm=='HEMOPTYSIS'~'COUGH',
                            saterm=='DIARRHEA'~'DIARRHOEA',
                            saterm=='CONJUNCTIVAL CONGESTION '~'CONJUNCTIVITIS',
                            saterm%like%'FEVER'~'HISTORY OF FEVER',
                            saterm=='SEIZURE'~'SEIZURES',
                            saterm%like%'TRANSPLANT'~'TRANSPLANTATION',
                            TRUE~saterm))%>%
    
    mutate(saterm=case_when(saterm%like%'ANOSMIA'~'LOST/ALTERED SENSE OF SMELL',
                            saterm%like%'AGEUSIA'~'LOST/ALTERED SENSE OF TASTE',
                            saterm=='LOSS OF SMELL (ANOSMIA)'~'LOST/ALTERED SENSE OF SMELL',
                            saterm=='LOSS OF TASTE (AGEUSIA)'~'LOST/ALTERED SENSE OF TASTE',
                            saterm=='LOST OF SMELL'~'LOST/ALTERED SENSE OF SMELL',
                            saterm=='LOST/ALTERED SENSE OF TASTE'~'LOST/ALTERED SENSE OF TASTE',
                            saterm=='LOST/ALTERED SENSE OF SMELL'~'LOST/ALTERED SENSE OF SMELL',
                            saterm=='DISTURBANCE OR LOSS OF SMELL (ANOSMIA)'~'LOST/ALTERED SENSE OF SMELL',
                            saterm=='DISTURBANCE OR LOSS OF TASTE (AGEUSIA)'~'LOST/ALTERED SENSE OF TASTE',
                            saterm=='ANOSMIA (LOSS OF SMELL OR TASTE)'~'LOST/ALTERED SENSE OF SMELL',
                            saterm=='AGEUSIA (LOSS OF TASTE)'~'LOST/ALTERED SENSE OF TASTE',
                            saterm=='LOSS OF SMELL (ANOSMIA)'~'LOST/ALTERED SENSE OF SMELL',
                            saterm=='LOSS OF TASTE (AGEUSIA)'~'LOST/ALTERED SENSE OF TASTE',
                            saterm=='LOSS OR DISTURBED SMELL'~'LOST/ALTERED SENSE OF SMELL',
                            saterm=='LOSS OR DISTURBED TASTE'~'LOST/ALTERED SENSE OF TASTE',
                            saterm=="NAUSEA/VOMITING"~'VOMITING/NAUSEA',
                            saterm%like%'MYALGIA OR FATIGUE'~'MUSCLE ACHES/JOINT PAIN',
                            saterm%like%'JOINT PAIN'~'MUSCLE ACHES/JOINT PAIN',
                            saterm%like%'MUSCLE ACHES'~'MUSCLE ACHES/JOINT PAIN',
                            saterm%like%'OTHER SIGN OR SYMPTOM'~'OTHER SIGNS AND SYMPTOMS',
                            saterm=='LOWER CHEST WALL INDRAWING'~'SHORTNESS OF BREATH',
                            saterm%like%'DEHYDRATION'~'SEVERE DEHYDRATION',
                            saterm%like%'RASH'~'SKIN RASH',
                            saterm=='EARPAIN'~'EAR PAIN',
                            TRUE ~ saterm))%>%
    
    mutate(saterm=case_when(saterm=='ADENOCARCINOMA'~'MALIGNANT NEOPLASM',
                            saterm=='ADENO-CARCINOMA'~'MALIGNANT NEOPLASM',
                            saterm=='ADENO CA'~'MALIGNANT NEOPLASM',
                            saterm=='ADENOCARCINOME'~'MALIGNANT NEOPLASM',
                            saterm=='ALL'~'MALIGNANT NEOPLASM',
                            saterm=='AML'~'MALIGNANT NEOPLASM',
                            saterm=='ASTROCYTOMA'~'MALIGNANT NEOPLASM',
                            saterm=='BCC'~'MALIGNANT NEOPLASM',
                            saterm=='BLASTOMA'~'MALIGNANT NEOPLASM',
                            saterm=='CA'~'MALIGNANT NEOPLASM',
                            saterm=='CANCER'~'MALIGNANT NEOPLASM',
                            saterm=='CARCINOMA'~'MALIGNANT NEOPLASM',
                            saterm=='CARCINOID'~'MALIGNANT NEOPLASM',
                            saterm=='CHOLANGIOCARCINOMA'~'MALIGNANT NEOPLASM',
                            saterm=='CLL'~'MALIGNANT NEOPLASM',
                            saterm=='CML'~'MALIGNANT NEOPLASM',
                            saterm=='DLBCL'~'MALIGNANT NEOPLASM',
                            saterm=='EPENDYOMA'~'MALIGNANT NEOPLASM',
                            saterm=='EWING'~'MALIGNANT NEOPLASM',
                            TRUE ~ saterm))%>%
    
    mutate(saterm=case_when(saterm=='GASTROINTESTINAL STROMAL TUMOUR'~'MALIGNANT NEOPLASM',
                            saterm=='GMB [NOT ANTI-GBM]'~'MALIGNANT NEOPLASM',
                            saterm=='GLIOBLASTOMA'~'MALIGNANT NEOPLASM',
                            saterm=='GIST'~'MALIGNANT NEOPLASM',
                            saterm=='HCC'~'MALIGNANT NEOPLASM',
                            saterm=='HODGKIN'~'MALIGNANT NEOPLASM',
                            saterm=='KAPOSI'~'MALIGNANT NEOPLASM',
                            saterm=='LEUKAEMIA'~'MALIGNANT NEOPLASM',
                            saterm=='LEUKEMIA'~'MALIGNANT NEOPLASM',
                            saterm=='LYMPHOMA'~'MALIGNANT NEOPLASM',
                            saterm=='MALIGNANCY'~'MALIGNANT NEOPLASM',
                            saterm=='MALIGNANT'~'MALIGNANT NEOPLASM',
                            saterm=='MDS'~'MALIGNANT NEOPLASM',
                            saterm=='MELANOMA'~'MALIGNANT NEOPLASM',
                            saterm=='MESOTHELIOMA'~'MALIGNANT NEOPLASM',
                            saterm=='MET'~'MALIGNANT NEOPLASM',
                            saterm=='METASTASES'~'MALIGNANT NEOPLASM',
                            saterm=='METASTATIC'~'MALIGNANT NEOPLASM',
                            saterm=='METS'~'MALIGNANT NEOPLASM',
                            saterm=='MYCOSIS FUNGOIDES'~'MALIGNANT NEOPLASM',
                            saterm=='MYELODYSPLASIA'~'MALIGNANT NEOPLASM',
                            saterm=='MYELODYSPLASTIC SYNDROME'~'MALIGNANT NEOPLASM',
                            saterm=='MYELOMA'~'MALIGNANT NEOPLASM',
                            TRUE ~ saterm))%>%
    
    mutate(saterm=case_when(saterm=='NHL'~'MALIGNANT NEOPLASM',
                            saterm=='NEUROBLASTOMA'~'MALIGNANT NEOPLASM',
                            saterm=='NEUROENDOCRINE CANCER'~'MALIGNANT NEOPLASM',
                            saterm=='NEUROENDOCRINE NEOPLASM'~'MALIGNANT NEOPLASM',
                            saterm=='NEUROENDOCRINE TUMOUR'~'MALIGNANT NEOPLASM',
                            saterm=='NON-HODGKIN'~'MALIGNANT NEOPLASM',
                            saterm=='RETINOBLASTOMA'~'MALIGNANT NEOPLASM',
                            saterm=='SARCOMA'~'MALIGNANT NEOPLASM',
                            saterm=='SCC'~'MALIGNANT NEOPLASM',
                            saterm=='SEZARY'~'MALIGNANT NEOPLASM',
                            saterm=='TCC'~'MALIGNANT NEOPLASM',
                            saterm=='WILMS'~'MALIGNANT NEOPLASM',
                            saterm=='CANER'~'MALIGNANT NEOPLASM',
                            saterm=='CACNER'~'MALIGNANT NEOPLASM',
                            saterm=='CANVER'~'MALIGNANT NEOPLASM', 
                            saterm=='CNACRE'~'MALIGNANT NEOPLASM',
                            saterm=='leukemia'~'MALIGNANT NEOPLASM', 
                            saterm=='leukaemia'~'MALIGNANT NEOPLASM',
                            TRUE ~ saterm))%>%
    
    mutate(saterm = iconv(saterm, to ="ASCII//TRANSLIT") %>% tolower()) %>%
    mutate(saterm = str_remove_all(saterm, "\\s*\\([^)]*\\)")) %>%
    mutate(saterm = str_replace_all(saterm, " - ", "_")) %>%
    mutate(saterm = str_replace_all(saterm, "/| / ", "_")) %>%
    mutate(saterm = str_replace_all(saterm, " ", "_")) %>%
    arrange(desc(saoccur))%>%
    distinct(usubjid,saterm, .keep_all =T)
  
  if(dtplyr.step){
    return(out)
  } else {
    return(out %>% as_tibble())
  }
}

imp_sa<-import.symptom.and.comorbidity.data(sa, dtplyr.step = FALSE)
save(imp_sa, file = "imp_sa.rda")


#' Process data on comorbidities
#' @param input Either the path of the symptoms/comorbidities data file (CDISC format) or output of \code{import.symptom.and.comorbidity.data}
#' @param dtplyr.step Return the output as \code{dtplyr_step} to avoid unnecessary future calls to \code{as_tibble} or \code{as.data.table}
#' @import dplyr tibble stringr tidyfast
#' @importFrom data.table as.data.table
#' @importFrom glue glue
#' @return Formatted comorbidity data as a tibble or \code{dtplyr_step}
#' @export process.comorbidity.data
process.comorbidity.data <- function(input,  minimum=100, dtplyr.step = FALSE){
  
  comorbid <- imp_sa%>%
    filter(sacat=="MEDICAL HISTORY") %>%
    filter(!is.na(sacat))%>%
    filter(!is.na(saterm))%>%
    arrange(desc(saoccur))%>%
    group_by(saterm) %>% 
    arrange(desc(saoccur))%>%
    mutate(n = sum(!is.na(saoccur))) %>%
    filter(n >= eval(!!minimum))%>%
    ungroup()%>%
    mutate(saterm = paste0("comorbid_",saterm)) %>%
    as.data.table() %>%
    pivot_wider(id_cols = usubjid, names_from = saterm, values_from = saoccur) 
  
  if(dtplyr.step){
    return(comorbid %>% lazy_dt(immutable = FALSE))
  } else {
    return(comorbid %>% as_tibble())
  }
}

imp_comorb<-process.comorbidity.data(imp_sa, minimum=100, dtplyr.step = FALSE)
save(imp_comorb, file = "imp_comorb.rda")

#' Process data on symptoms
#' @param input Either the path of the symptoms/comorbidities data file (CDISC format) or output of \code{import.symptom.and.comorbidity.data}
#' @param dtplyr.step Return the output as \code{dtplyr_step} to avoid unnecessary future calls to \code{as_tibble} or \code{as.data.table}
#' @import dplyr tibble tidyfast dtplyr
#' @importFrom data.table as.data.table
#' @importFrom glue glue
#' @return Formatted symptom data as a tibble or \code{dtplyr_step}
#' @export process.symptom.data
process.symptom.data <- function(input,  minimum=100, dtplyr.step = FALSE){
  
  symptom_w <- imp_sa%>% mutate(studyid=substr(usubjid,1, 7))%>%
    filter(studyid!="CVZXZMV")%>%
    filter(sacat=="SIGNS AND SYMPTOMS AT HOSPITAL ADMISSION" | 
             (sacat=="SIGNS AND SYMPTOMS AT INITIAL ACUTE COVID-19 ILLNESS" & fup_hospitalised=="TRUE"))%>% #Add followup patients with acute hospitalisation
    filter(saterm!="covid-19_symptoms")%>%
    arrange(desc(saoccur))%>%
    group_by(saterm) %>% 
    arrange(desc(saoccur))%>%
    mutate(n = sum(!is.na(saoccur))) %>%
    filter(n >= eval(!!minimum))%>%
    ungroup()%>%
    mutate(saterm = paste0("symptoms_",saterm)) %>%
    as.data.table() %>%
    pivot_wider(id_cols = usubjid, names_from = saterm, values_from = saoccur) %>%
    as.data.frame()
  
  date_onset<-imp_sa%>%
    ungroup()%>%
    filter(sacat=="SIGNS AND SYMPTOMS AT HOSPITAL ADMISSION")%>%
    filter(saoccur==TRUE) %>%
    distinct(usubjid, .keep_all =T)%>%
    select(usubjid, "date_onset"=sastdy)
  
  covid_clinic_diagn<- imp_sa%>%
    filter(sacat=="CLINICALLY-DIAGNOSED COVID-19")%>%
    mutate(saoccur=case_when(is.na(sapresp)~TRUE,
                             TRUE~saoccur))%>%
    arrange(desc(saoccur))%>%
    distinct(usubjid, .keep_all =T)%>%
    rename("clin_diag_covid_19"=saoccur)%>%
    select(usubjid,clin_diag_covid_19)
  
  symptomatic<-imp_sa%>%mutate(studyid=substr(usubjid,1, 7))%>%
    filter(studyid!="CVZXZMV")%>%
    ungroup()%>%
    filter(sacat=="SIGNS AND SYMPTOMS AT HOSPITAL ADMISSION")%>%
    mutate(symptomatic=case_when(saterm=="asymptomatic" & saoccur==TRUE~FALSE,
                                 saterm=="asymptomatic" & saoccur==FALSE~TRUE,
                                 TRUE~saoccur))%>%
    arrange(desc(symptomatic))%>%
    distinct(usubjid, .keep_all =T)%>%
    select(usubjid, symptomatic)
  
  symptom<- date_onset%>%
    full_join(covid_clinic_diagn, by=c("usubjid"))%>%
    full_join(symptomatic, by = c("usubjid"))%>%
    full_join(symptom_w, by = c("usubjid"))
  
  if(dtplyr.step){
    return(symptom %>% lazy_dt(immutable = FALSE))
  } else {
    return(symptom %>% as_tibble())
  }
}

imp_symptom<-process.symptom.data(imp_sa, minimum=100, dtplyr.step = FALSE)
save(imp_symptom, file = "imp_symptom.rda")


#' Process data on treatments
#' @param file.name Path of the intervention data file (CDISC format)
#' @param dtplyr.step Return the output as \code{dtplyr_step} to avoid unnecessary future calls to \code{as_tibble} or \code{as.data.table}
#' @import dplyr tibble stringr
#' @return Formatted treatment data (long format) as a tibble or \code{dtplyr_step}
#' @export process.treatment.data

#Import IN data file
int <- fread(file=" Partner_IN_20230110", skip = 0) 
int<-as_tibble(int)
colnames(int) <- tolower (colnames(int))

process.treatment.data <- function(file.name,  dtplyr.step = FALSE){
  
  treatment<-int%>%
    filter(inpresp =="Y") %>%
    filter(inevintx!="BEFORE HOSPITAL ADMISSION")%>%
    mutate(inoccur = case_when(inoccur == "Y" ~ TRUE,
                               inoccur == "N" ~ FALSE,
                               TRUE ~ NA))%>%
    filter(!is.na(inoccur))%>%
    filter(incat!="MEDICAL HISTORY" | is.na (incat))%>%
    mutate(intrt_original=intrt)%>%
    mutate(intrt=toupper(intrt))%>%
    mutate(intrt=as.character(intrt))%>%
    mutate(inmodify=as.character(inmodify))%>%
    mutate(incat=as.character(incat))%>%
    mutate(intrt=case_when(inmodify!=""~inmodify,
                           TRUE ~ intrt))%>%
    mutate(intrt=case_when(incat=="EXTRACORPOREAL"~'EXTRACORPOREAL',
                           incat=="INVASIVE VENTILATION"~'INVASIVE VENTILATION',
                           incat=="INVASIVE VENTILATION"~'INVASIVE VENTILATION',
                           incat=="NON-INVASIVE VENTILATION "~'NON-INVASIVE VENTILATION ',
                           incat=="OTHER INTEVENTIONS"~'OTHER INTERVENTIONS',
                           incat=="ANTIBIOTIC AGENTS"~ "ANTIBIOTIC AGENTS",
                           incat=="ANTIFUNGAL AGENTS"~ "ANTIFUNGAL AGENTS",
                           incat=="ANTIVIRAL AGENTS"~ "ANTIVIRAL AGENTS",
                           incat=="CORTICOSTEROIDS"~ "CORTICOSTEROIDS",
                           incat=="ANTIMALARIAL AGENTS"~ "ANTIMALARIAL AGENTS",
                           incat=="NSAIDS"~"NON-STEROIDAL ANTI-INFLAMMATORY (NSAIDS)",
                           TRUE~intrt)) %>%
    mutate(intrt=case_when(intrt%like%'ECMO'~'EXTRACORPOREAL',
                           intrt=='EXTRA CORPOREAL LIFE SUPPORT'~'EXTRACORPOREAL',
                           intrt=='EXTRACORPOREAL SUPPORT'~'EXTRACORPOREAL',
                           intrt=='CONTINUOUS RENAL REPLACEMENT THERAPIES (CRRT)'~'RENAL REPLACEMENT THERAPIES',
                           intrt%like%'RENAL REPLACEMENT THERAPY' |
                             intrt%like% 'DIALYSIS'~ 'RENAL REPLACEMENT THERAPIES',
                           intrt%like% 'HEMOFILTRATION'~ 'RENAL REPLACEMENT THERAPIES',
                           intrt=='ERP CVVH'~ 'RENAL REPLACEMENT THERAPIES',
                           TRUE ~ intrt))%>%
    ###IMV
    mutate(intrt=case_when(intrt=='INVASIVE MECHANICAL LUNG VENTILATION'~'INVASIVE VENTILATION',
                           intrt=='INVASIVE MECHANICAL VENTILATION'~'INVASIVE VENTILATION',
                           intrt=='MECHANICAL VENTILATION'~'INVASIVE VENTILATION',
                           intrt=='RE-INTUBATION'~'INVASIVE VENTILATION',
                           intrt=='INVASIVE VENTILATION'~'INVASIVE VENTILATION',
                           intrt%like%'APRV'~'INVASIVE VENTILATION',
                           intrt=='INTUBATION AND MECHANICAL VENTILATION'~'INVASIVE VENTILATION',
                           intrt=='MECHANICAL SUPPORT'~'INVASIVE VENTILATION',
                           intrt%like%'EXTUBATION'~'INVASIVE VENTILATION',
                           intrt=="VENTILATED"~'INVASIVE VENTILATION',
                           TRUE ~ intrt))%>%
    
    ###NIV
    mutate(intrt=case_when(intrt%like%'CPAP'~'NON-INVASIVE VENTILATION',
                           intrt%like%'BIPAP'~'NON-INVASIVE VENTILATION',
                           intrt%like%'NON-INVASIVE MECHANICAL VENTILATION (BIPAP, CPAP, OCNAF (OPTIFLOW) ...)'~'NON-INVASIVE VENTILATION',
                           intrt%like%'NON-INVASIVE VENTILATION'~'NON-INVASIVE VENTILATION',
                           intrt=='NON-INVASIVE MECHANICAL VENTILATION'~'NON-INVASIVE VENTILATION',
                           intrt=='NON-INVASIVE POSITIVE PRESSURE VENTILATION'~'NON-INVASIVE VENTILATION',
                           intrt=='NON-INVASIVE RESPIRATORY SUPPORT'~'NON-INVASIVE VENTILATION',
                           TRUE ~ intrt))%>%
    
    mutate(intrt=case_when(intrt%like%'OTHER INTERVENTION'~'OTHER INTERVENTIONS',
                           intrt%like%'CHEMOTHERAPY'| intrt%like%'ANTI-DIABETIC MEDICATIONS'|intrt%like%'BRONCHOSCOPY'|
                             intrt%like%'PROTON PUMP INHIBITORS'|intrt%like%'STATINS'|intrt%like%'MORPHINE'|
                             intrt%like%'HALOPERIDOL'|intrt%like%'OLANZAPINE'~'OTHER INTERVENTIONS',
                           intrt=='OTHER TARGETED COVID-19 MEDICATIONS'~'OTHER INTERVENTIONS',
                           intrt=='OTHER TREATMENTS FOR COVID19'~'OTHER INTERVENTIONS',
                           intrt%like%"NON-STEROIDAL"~"NON-STEROIDAL ANTI-INFLAMMATORY",
                           intrt%like%"NON STEROIDAL"~"NON-STEROIDAL ANTI-INFLAMMATORY",                           
                           TRUE ~ intrt))%>%
    
    mutate(intrt=case_when(intrt%like%'SURGICAL FEEDING TUBE'~'TOTAL PARENTERAL NUTRITION',
                           
                           #MASKS 
                           intrt=="FACE MASK"~'MASK OXYGEN THERAPY',
                           intrt=="NASAL CANNULA"~'NASAL OXYGEN THERAPY',
                           
                           ####HFNC 
                           intrt=='OXYGEN THERAPY WITH HIGH FLOW NASAL CANULA'~'HIGH-FLOW NASAL CANULA OXYGEN THERAPY',
                           intrt=='HIGH-FLOW NASAL CANNULA OXYGEN THERAPY'~'HIGH-FLOW NASAL CANULA OXYGEN THERAPY',
                           intrt=='HUMIDIFIED HIGH FLOW NASAL CANNULA (HHFNC)'~'HIGH-FLOW NASAL CANULA OXYGEN THERAPY',
                           intrt=='HIGH-FLOW NASAL OXYGEN'~'HIGH-FLOW NASAL CANULA OXYGEN THERAPY',
                           intrt=='AIRVO (HIGH FLOW NASAL CANULA)'~'HIGH-FLOW NASAL CANULA OXYGEN THERAPY',
                           
                           ###PRONE POSITIONING 
                           intrt%like%'PRONACI'~'PRONE POSITIONING',
                           intrt=='PRONE POSITIONING'~'PRONE POSITIONING',
                           
                           ###TRACHEOSTOMY 
                           intrt%like%'TRACHEOSTOMY'~'TRACHEOSTOMY',
                           
                           ###NITRIC OXIDE 
                           intrt%like%'NITRIC OXIDE'~'INHALED NITRIC OXIDE',
                           intrt%like%"NITROUS OXIDE" ~ "INHALED NITRIC OXIDE",
                           
                           ###RESPIRATORY SUPPORT 
                           intrt=='RESPIRATORY SUPPORT'~'RESPIRATORY SUPPORT',
                           
                           ###Corticosteroids
                           intrt=="CORTICOSTEROID"~ "CORTICOSTEROIDS",
                           intrt=="DEXAMETHASONE"~ "CORTICOSTEROIDS",
                           intrt=="BETAMETHASONE"~ "CORTICOSTEROIDS",
                           intrt%like%"PREDNISOLONE"~ "CORTICOSTEROIDS",
                           intrt=="ORAL STEROIDS"~ "CORTICOSTEROIDS",
                           intrt=="STEROIDS"~ "CORTICOSTEROIDS",
                           intrt%like%"HYDROCORTISONE"~ "CORTICOSTEROIDS",
                           intrt%like%"BLOOD TRANSFUSION OR BLOOD PRODUCT"~ "BLOOD TRANSFUSION OR BLOOD PRODUCT",
                           TRUE ~ intrt))%>%
    
    mutate(intrt=case_when(intrt%like%"ANTIVIRAL" ~ "ANTIVIRAL AGENTS",
                           intrt%like%"ARV" ~ "ANTIVIRAL AGENTS",
                           intrt%like%"ANTIRETROVIRAL" ~ "ANTIVIRAL AGENTS",
                           intrt%like%"RIBAVIRIN" ~ "ANTIVIRAL AGENTS",
                           intrt%like%"LOPINAVIR AND RITONAVIR" ~ "ANTIVIRAL AGENTS",
                           intrt%like%"LOPINAVIR" ~ "ANTIVIRAL AGENTS",
                           intrt%like%"OSELTAMIVIR" ~ "ANTIVIRAL AGENTS",
                           intrt%like%"REMDESIVIR" ~ "ANTIVIRAL AGENTS",
                           intrt%like%"NEURAMINIDASE INHIBITORS" ~ "ANTIVIRAL AGENTS",
                           intrt%like%"ZANAMIVIR" ~ "ANTIVIRAL AGENTS",
                           intrt%like%"RIBAVARIN" ~ "ANTIVIRAL AGENTS",
                           TRUE ~ intrt))%>%
    
    #mutate(intrt=case_when(intrt%like%"FAVIPIRAVIR" ~ "ANTIVIRAL AGENTS",
    #                      intrt%like%"ATAZANAVIR" ~ "ANTIVIRAL AGENTS",
    #                     TRUE ~ intrt))%>%
    
    mutate(intrt=case_when(intrt%like%"FLUCLOXACILLIN"~ "ANTIBIOTIC AGENTS",
                           intrt%like%"ANTIBIOTIC"~ "ANTIBIOTIC AGENTS",
                           intrt%like%"AMIKACIN"~ "ANTIBIOTIC AGENTS",
                           intrt%like%"AMOX"~ "ANTIBIOTIC AGENTS",
                           intrt%like%"AUGUMENTIN"~ "ANTIBIOTIC AGENTS",
                           intrt%like%"AZITHROMYCIN"~ "ANTIBIOTIC AGENTS",
                           intrt%like%"AZITHRYOMYCIN"~ "ANTIBIOTIC AGENTS",
                           intrt%like%"BENZY"~ "ANTIBIOTIC AGENTS",
                           intrt%like%"AUGUMENTIN"~ "ANTIBIOTIC AGENTS",
                           TRUE ~ intrt))%>%
    
    mutate(intrt=case_when(intrt%like%"AZITHRYOMYCIN"~ "ANTIBIOTIC AGENTS",
                           intrt%like%"CEFTR"~ "ANTIBIOTIC AGENTS",
                           intrt%like%"CEFR"~ "ANTIBIOTIC AGENTS",
                           intrt%like%"DOXYCYCLINE"~ "ANTIBIOTIC AGENTS",
                           intrt%like%"CHLORAMPHENICOL"~ "ANTIBIOTIC AGENTS",
                           intrt%like%"CIPROFLOXACIN"~ "ANTIBIOTIC AGENTS",
                           intrt%like%"GENTAMICIN"~ "ANTIBIOTIC AGENTS",
                           intrt%like%"MEROPENEM"~ "ANTIBIOTIC AGENTS",
                           intrt%like%"METRONIDAZOLE"~ "ANTIBIOTIC AGENTS",
                           TRUE ~ intrt))%>%
    
    mutate(intrt=case_when(intrt%like%"ANTIMALARIAL" | intrt%like%"CHLOROQUINE" ~ "ANTIMALARIAL AGENTS",
                           intrt%like%"ANTIFUNGAL" ~ "ANTIFUNGAL AGENTS",
                           intrt %like% "OROGASTRIC"~"NASO/ NASOGASTRIC ORAL/OROGASTRIC FLUIDS",
                           intrt %like% "NGT OR OGT REQUIRED FOR NUTRITION"~"NASO/ NASOGASTRIC ORAL/OROGASTRIC FLUIDS",
                           intrt%like%'DOBUTAMINE' |  intrt%like%'DOPAMINE' |  intrt%like%'MILRINONE' 
                           |  intrt%like%'LEVOSIMENDAN' |  intrt%like%'EPINEPHRINE' |  intrt%like%'NOREPINEPRINE'
                           |  intrt%like%'INOTROPES' |intrt%like%'VASOPRESS' |intrt%like%'NORADRENALINE' |
                             intrt%like%'ADRENALINE' |intrt%like%'BETA BLOCKER' ~'INOTROPES / VASOPRESSORS',
                           TRUE ~ intrt))%>%
    
    mutate(intrt=case_when(intrt%like%'IMMUNOGLOBULI' ~ "IMMUNOGLOBULI",
                           intrt=='CONVALESCENT PLASMA' ~ "CONVALESCENT PLASMA",
                           intrt%like%'IMMUNOSUPPRES' ~ "IMMUNOSUPPRESSANTS",
                           intrt%like%'IMMUNOSTIMULANTS' ~ "IMMUNOSUPPRESSANTS",
                           intrt%like%'IMMUNOTHERAPY' ~ "IMMUNOSUPPRESSANTS",
                           intrt=="IL6 INHIBITOR" ~ "IMMUNOSUPPRESSANTS",
                           intrt=="TOCILIZUMAB" ~ "IMMUNOSUPPRESSANTS",
                           intrt%like%"INTERFERON" ~ "IMMUNOSTIMULANTS",
                           intrt%like%"HEPARIN" ~ "THERAPEUTIC ANTICOAGULANT",
                           intrt%like%"NOXAPARIN" ~ "THERAPEUTIC ANTICOAGULANT",
                           intrt=="ENOXAPARIN" ~ "THERAPEUTIC ANTICOAGULANT",
                           TRUE ~ intrt))%>%
    
    mutate(intrt=case_when(intrt%like%"SPIRONOLACTONE" ~ "DIURETICS",
                           intrt%like%"DIURETIC" ~ "DIURETICS",
                           #intrt%like%"NITROUS OXIDE" ~ "inhaled_nitric_oxide",
                           intrt=="CPR" ~ "Cardiopulmonary resuscitation",
                           intrt%like%"EXPERIMENTAL AGENT" ~ "EXPERIMENTAL AGENTS",
                           intrt%like%"SARILUMAB" ~ "EXPERIMENTAL AGENTS",
                           intrt%like%"IV FLUID" ~ "INTRAVENOUS FLUIDS",
                           intrt%like%"I.V. SOLUTIONS" ~ "INTRAVENOUS FLUIDS",
                           intrt %like% "ANGIOTENSIN" | intrt %like% "ACE"~ "AGENTS ACTING ON THE RENIN-ANGIOTENSIN SYSTEM",
                           intrt%like%"ANTIINFLAMMATORY" ~ "ANTIINFLAMMATORY",
                           TRUE ~ intrt))%>%    as.data.frame()%>%
    select(studyid,usubjid,'treatment'=intrt,inoccur,intrt_original,inmodify,incat, inevintx, indur,instdy,incdstdy, inendy,indy,instdy,inindc)%>%
    mutate(treatment = iconv(treatment, to ="ASCII//TRANSLIT") %>% tolower()) %>%
    mutate(treatment = str_remove_all(treatment, "\\s*\\([^)]*\\)")) %>%
    mutate(treatment = str_replace_all(treatment, " - ", "_")) %>%
    mutate(treatment = str_replace_all(treatment, "-", "_")) %>%
    mutate(treatment = str_replace_all(treatment, "/| / ", "_")) %>%
    mutate(treatment = str_replace_all(treatment, " ", "_"))
  
  if(dtplyr.step){
    return(treatment)
  } else {
    return(treatment %>% as_tibble())
  }
}

imp_int<-process.treatment.data(int, dtplyr.step = FALSE)
save(imp_int, file = "imp_int.rda")

#' Process data on the most common treatments
#' @param input Either the path of the interventions data file (CDISC format) or output of \code{process.treatment.data}
#' @param minimum The minimum number of times a treatment need appear to be considered "common"; default 1000.
#' @param dtplyr.step Return the output as \code{dtplyr_step} to avoid unnecessary future calls to \code{as_tibble} or \code{as.data.table}
#' @import dplyr tibble dtplyr tidyfast
#' @importFrom data.table as.data.table
#' @importFrom glue glue
#' @return Formatted common treatment data (wide format) as a tibble or \code{dtplyr_step}
#' @export process.common.treatment.data
process.common.treatment.data <- function(file.name, minimum=10, dtplyr.step = FALSE){
  
  oxy_within_d1<-imp_int%>%
    mutate(treatment=case_when(treatment=="extracorporeal" | 
                                 treatment=="respiratory_support" |
                                 treatment=="high_flow_nasal_cannula" |
                                 treatment=="invasive_ventilation" |
                                 treatment=="mask_oxygen_therapy" |
                                 treatment=="nasal_oxygen_therapy" |
                                 treatment=="oxygen_therapy"|
                                 treatment=="non_invasive_ventilation"~"oxygen_therapy",
                               TRUE~treatment))%>%
    filter(treatment=="oxygen_therapy")%>%
    filter(inevintx=="AT HOSPITAL ADMISSION"|indy==1)%>%
    mutate(oxytreat_when=case_when(inevintx=="AT HOSPITAL ADMISSION"~"at_admi",
                                   indy==1~"within_24h"))%>%
    arrange(desc(inoccur))%>%
    distinct(usubjid,treatment, .keep_all =T)%>%
    select(usubjid,"d1_oxygen_therapy"=inoccur)
  
  treatment <- imp_int%>%
    group_by(treatment)%>% 
    arrange(desc(inoccur))%>%
    mutate(n = sum(!is.na(inoccur)))%>%
    filter(n >= eval(!!minimum))%>%
    ungroup()%>%
    filter(treatment!="extracorporeal" & 
             treatment!="inhaled_nitric_oxide" &
             treatment!="oxygen_therapy" &
             treatment!="prone_position_ventilation" &
             treatment!="prone_ventilation" &
             treatment!="respiratory_support" &
             treatment!="tracheostomy" &
             treatment!="prone_positioning")%>%
    filter(treatment!="covid_19_vaccination")%>%
    filter(treatment!="supplemental_oxygen_fio2")%>%
    arrange(desc(inoccur))%>%
    distinct(usubjid, treatment, .keep_all =T)%>% 
    mutate(treatment = glue("treat_{treatment}", treatment = treatment))%>%
    as.data.table()%>%
    dt_pivot_wider(id_cols = usubjid, names_from = treatment,  values_from = inoccur)%>%
    as.data.frame()
  
  ####calculating oxygen therapy overall
  treat_oxy <- imp_int%>%
    mutate(treatment=case_when(treatment=="extracorporeal" | 
                                 treatment=="respiratory_support" |
                                 treatment=="high_flow_nasal_cannula" |
                                 treatment=="invasive_ventilation" |
                                 treatment=="mask_oxygen_therapy" |
                                 treatment=="nasal_oxygen_therapy" |
                                 treatment=="oxygen_therapy"|
                                 treatment=="non_invasive_ventilation"~"treat_oxygen_therapy",
                               TRUE~treatment))%>%
    filter(treatment=="treat_oxygen_therapy")%>%
    arrange(desc(inoccur))%>%
    distinct(usubjid, .keep_all =T)%>%
    select(usubjid,"treat_oxygen_therapy"=inoccur)
  
  ###adding duration for invasive_ventilation and non_invasive_ventilation
  indur <- imp_int%>%
    select(usubjid,treatment, inoccur,indur,indy)%>%
    filter(treatment=="invasive_ventilation"|treatment=="non_invasive_ventilation")%>%
    select(usubjid,treatment, inoccur,indur,indy)%>%
    filter(treatment=="invasive_ventilation"|treatment=="non_invasive_ventilation")%>%
    mutate(treatment=case_when(treatment=='non_invasive_ventilation'~'dur_niv',
                               treatment=='invasive_ventilation'~'dur_imv',
                               TRUE~treatment))%>%
    mutate(indur_clean=as.numeric(gsub("[^0-9.]", "",indur)))%>%
    filter(!is.na(indur_clean)  | indur_clean!="")%>%
    distinct(usubjid,treatment, .keep_all =T)%>%
    pivot_wider(id_cols = usubjid, names_from = treatment,  values_from = indur_clean)%>%
    as_tibble()
  
  treatment <-treatment%>%
    full_join(treat_oxy)%>%
    full_join(oxy_within_d1)%>%
    full_join(indur)
  
  if(dtplyr.step){
    return(treatment) %>% lazy_dt(immutable = FALSE)
  } else {
    return(treatment %>% as_tibble())
  }
  
}

imp_treat<-process.common.treatment.data(imp_int, minimum=10, dtplyr.step = FALSE)
save(imp_treat, file = "imp_treat.rda")

#' Process data on the most common icu treatments
#' @param input Either the path of the interventions data file (CDISC format) or output of \code{process.treatment.data}
#' @param minimum The minimum number of times a treatment need appear to be considered "common"; default 1000.
#' @param dtplyr.step Return the output as \code{dtplyr_step} to avoid unnecessary future calls to \code{as_tibble} or \code{as.data.table}
#' @import dplyr tibble dtplyr tidyfast
#' @importFrom data.table as.data.table
#' @importFrom glue glue
#' @return Formatted common treatment data (wide format) as a tibble or \code{dtplyr_step}
#' @export process.treatment.icu.data
process.treatment.icu.data <- function(file.name,imp_icu,imp_dm,imp_ds, minimum=10, dtplyr.step = FALSE){
  
  adm_date<-imp_dm%>%
    select(usubjid,date_admit)
  
  icu_ever<-imp_icu%>%
    filter(ever_icu==TRUE)%>%
    left_join(adm_date)%>%
    left_join(imp_ds)%>%
    select(usubjid,hostdy,hoendy)
  
  treat_oxy_icu <- imp_int%>%
    mutate(treatment=case_when(treatment=="extracorporeal" | 
                                 treatment=="respiratory_support" |
                                 treatment=="high_flow_nasal_cannula" |
                                 treatment=="invasive_ventilation" |
                                 treatment=="mask_oxygen_therapy" |
                                 treatment=="nasal_oxygen_therapy" |
                                 treatment=="oxygen_therapy"|
                                 treatment=="non_invasive_ventilation"~"treat_oxygen_therapy",
                               TRUE~treatment))%>%
    
    filter(treatment=="treat_oxygen_therapy")%>%
    arrange(desc(inoccur))%>%
    left_join(icu_ever,by = c("usubjid"))%>%
    mutate(indy=as.numeric(indy))%>%
    mutate(hostdy=as.numeric(hostdy))%>%
    mutate(hoendy=as.numeric(hoendy))%>%
    mutate(int_icu=case_when((indy>=hostdy)~ TRUE, 
                             TRUE ~ FALSE))%>%
    filter(int_icu==TRUE)%>%
    arrange(desc(inoccur))%>%
    distinct(usubjid, treatment, .keep_all =T)%>%
    select(usubjid,"icu_treat_oxygen_therapy"=inoccur)
  
  imp_treat_icu<-imp_int%>%
    filter(!is.na(indy))%>%
    group_by(treatment) %>% 
    arrange(desc(inoccur))%>%
    mutate(n = sum(!is.na(inoccur))) %>%
    filter(n >= eval(!!minimum)) %>%
    ungroup()%>%
    filter(treatment!="covid_19_vaccination")%>%
    filter(treatment!="supplemental_oxygen_fio2")%>%
    filter(treatment!="extracorporeal" & 
             treatment!="inhaled_nitric_oxide" &
             treatment!="oxygen_therapy" &
             treatment!="prone_position_ventilation" &
             treatment!="prone_ventilation" &
             treatment!="respiratory_support" &
             treatment!="tracheostomy" &
             treatment!="prone_positioning")%>%
    arrange(desc(inoccur))%>%
    left_join(icu_ever,by = c("usubjid"))%>%
    mutate(indy=as.numeric(indy))%>%
    mutate(hostdy=as.numeric(hostdy))%>%
    mutate(hoendy=as.numeric(hoendy))%>%
    mutate(int_icu=case_when((indy>=hostdy)~ TRUE, 
                             TRUE ~ FALSE))%>%
    filter(int_icu==TRUE)%>%
    arrange(desc(inoccur))%>%
    distinct(usubjid, treatment, .keep_all =T)%>%
    mutate(treatment = glue("icu_treat_{treatment}", treatment = treatment)) %>%
    as.data.table() %>%
    pivot_wider(id_cols = usubjid, names_from = treatment,  values_from = inoccur)%>%
    full_join(treat_oxy_icu)
  
  if(dtplyr.step){
    return(imp_treat_icu) %>% lazy_dt(immutable = FALSE)
  } else {
    return(imp_treat_icu %>% as_tibble())
  }
  
}

imp_treat_icu<-process.treatment.icu.data(imp_int, imp_icu,imp_dm,imp_ds, minimum=10,dtplyr.step = FALSE)
save(imp_treat_icu, file = "imp_treat_icu.rda")


#Join all data =================================================================
import.tbl_feb2023.v1<-imp_dm%>%
  left_join(imp_mb, by = c("usubjid"))%>%
  left_join(imp_comorb, by = c("usubjid"))%>%
  left_join(imp_rp, by = c("usubjid"))%>%
  left_join(imp_symptom, by = c("usubjid"))%>%
  left_join(imp_treat, by = c("usubjid"))%>%
  left_join(imp_icu, by = c("usubjid"))%>%
  left_join(imp_treat_icu, by = c("usubjid"))%>%
  left_join(imp_lb, by = c("usubjid"))%>%
  left_join(imp_vs, by = c("usubjid"))%>%
  left_join(imp_ds, by = c("usubjid"))
save(import.tbl_feb2023.v1, file = "import.tbl_feb2023.v1.rda")

