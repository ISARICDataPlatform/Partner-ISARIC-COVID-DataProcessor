
# =========== THIS IS THE 2ND PART OF THE CODE. THIS IS RUN AFTER RUNNING THE CODE LABELLED "Partner_ImportFunctionsModular.R"=====
# =================================================================================================================================
  
  
#' Preprocessing step for all aggregations. Currently: remaps outcome to death, discharge or NA, cuts age into 5-year age groups, and adds a year-epiweek column
#' @param input.tbl Input tibble (output of \code{process.all.data})
#' @import dtplyr dplyr purrr lubridate tibble
#' @importFrom glue glue
#' @return A \code{tibble} intended for input into other aggregation functions (e.g. \code{age.pyramid.prep})
#' @export data.preprocessing


input.tbl<- import.tbl_partner

### PART 1 - Run from line 17 to 162
data.preprocessing <- function(input.tbl){
    
  #create upper respiratory tract symptoms combining several symptoms
  mutate(symptrcd_upper_respiratory_tract_symptoms=NA)%>%
  mutate(symptrcd_upper_respiratory_tract_symptoms=case_when(
    symptoms_upper_respiratory_tract_symptoms==FALSE|
      symptoms_sore_throat==FALSE|
      symptoms_runny_nose==FALSE|
      symptoms_ear_pain==FALSE~FALSE,
    TRUE~symptrcd_upper_respiratory_tract_symptoms))%>%
  mutate(symptrcd_upper_respiratory_tract_symptoms=case_when(
    symptoms_upper_respiratory_tract_symptoms==TRUE|
      symptoms_sore_throat==TRUE|
      symptoms_runny_nose==TRUE|
      symptoms_ear_pain==TRUE~TRUE,
    TRUE~symptrcd_upper_respiratory_tract_symptoms))%>%
  mutate(age=replace(age,age>120,NA))%>%
  mutate(agegp10 = cut(age, right = FALSE, breaks = c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 120))) %>%
  mutate(agegp5 = cut(age, right = FALSE, breaks = c(0,5, 10,15, 20,25, 30,35, 40,45, 50,55,
                                                     60,65, 70,75, 80,85, 90, 95, 100, 120))) %>%

    mutate(calendar.year.admit = year(date_admit)) %>%
    mutate(calendar.month.admit = month(date_admit)) %>%
    mutate(slider_monthyear = map2_chr(calendar.year.admit, calendar.month.admit, month.year.mapper)) %>%
    mutate(year.admit = map_dbl(date_admit, epiweek.year)) %>%
    mutate(epiweek.admit = epiweek(date_admit)) %>%
    mutate(year.epiweek.admit=paste0(year.admit,"-", epiweek.admit))%>%
    mutate(year.epiweek.admit = replace(year.epiweek.admit, year.epiweek.admit == "NA-NA", NA)) %>%
    mutate(lower.age.bound  = map_dbl(agegp10, extract.age.boundaries, TRUE)) %>%
    mutate(upper.age.bound  = map_dbl(agegp10, extract.age.boundaries, FALSE)) %>%
    mutate(slider_agegp10 = fct_relabel(agegp10, prettify.age.labels)) %>%
    dplyr::select(-agegp10) %>%
    
    #rename slider variables
    rename(slider_icu_ever = ever_icu) %>%
    rename(slider_country = country) %>%
    rename(slider_sex = sex) %>%
    rename(slider_symptomatic = symptomatic) %>%
    
    #set as NA implausible negative value 
    mutate_at(vars(all_of(c(starts_with("t_"),starts_with("dur_")))), function(x){replace(x,x<0,NA)})%>%
    
    #deleting implausible respiratory rates based on age
    mutate(vs_resp=case_when(vs_resp<= 3 ~ NA_real_,
                             vs_resp<=5 & age < 10 ~ NA_real_ ,
                             TRUE ~ vs_resp)) %>% 
    
    ##################################################################
  #set as NA outliers for vital sign and laboratory variables
  mutate(age_outlier = ifelse(age>10,1,0))%>% 
    group_by(age_outlier)%>% 
    mutate_at(vars(c(all_of(c(starts_with("vs_"),starts_with("lab_"))))), 
              function(x,na.rm = FALSE){replace(x, 
                                                x<(quantile(x, 0.025, na.rm = TRUE))|
                                                  x>(quantile(x, 0.975, na.rm = TRUE)),
                                                NA_real_)
              })%>% 
    ungroup()%>%
    mutate_at(vars(c(all_of(c(starts_with("t_"),starts_with("dur_"))))), 
              function(x,na.rm = FALSE){replace(x, 
                                                x>(quantile(x, 0.975, na.rm = TRUE)),
                                                NA_real_)
              })%>% 
  
  #############################################################################################
  #Calculating BMI
  
  mutate(vs_bmi_calc=vs_weight/(vs_height/100)^2)%>%
    mutate(vs_bmi_calc=as.numeric(vs_bmi_calc))%>%
    mutate(vs_bmi=as.numeric(vs_bmi))%>%
    mutate(bmi_comb=ifelse(!is.na(vs_bmi),vs_bmi,vs_bmi_calc))%>%
    mutate(und_nutr=case_when(bmi_comb<18.5 & age<65~"under nutrition",
                              bmi_comb<20.5 & age>64~"under nutrition",
                              bmi_comb>18.4 & age<65 ~"normal nutrition",
                              bmi_comb>20.4 & age>65 ~"normal nutrition",
                              TRUE~NA_character_))%>%
    mutate(embargo_length=case_when(date_admit>date_pull-14~TRUE,
                                    date_admit<=date_pull-14~FALSE
    ))%>%
    
    as_tibble()
}


#' @keywords internal
#' @export prettify.age.labels
prettify.age.labels <- function(a){
  temp <- substr(a, 2, nchar(a) - 1)
  newlabels <- map_chr(temp, function(x) {
    components <- as.numeric(str_split_fixed(x, ",", Inf))
    components[2] <- components[2] - 1
    paste(components, collapse = "-")
  })
  str_replace(newlabels, "90-119", "90+")
}

#' @keywords internal
#' @export extract.age.boundaries
extract.age.boundaries <- function(agestring, lower = TRUE){
  agestring <- as.character(agestring)
  temp <- substr(agestring, 2, nchar(agestring)-1)
  if(lower){
    as.numeric(str_split_fixed(temp, ",", Inf)[1])
  } else {
    as.numeric(str_split_fixed(temp, ",", Inf)[2]) - 1
  }
}

#' @keywords internal
#' @export cleaning.unplosible.dates
cleaning.unplosible.dates <- function(date){
  if(is.na(date)){
    return(NA)
  }
  if(year(date)==2019 & date > ymd("2019-12-28")){
    2020
  } else {
    year(date)
  }
}


#' @keywords internal
#' @export epiweek.year
epiweek.year <- function(date){
  if(is.na(date)){
    return(NA)
  }
  if(year(date)==2019 & date > ymd("2019-12-28")){
    2020
  } else {
    year(date)
  }
}

#' @keywords internal
#' @export month.year.mapper
month.year.mapper <- function(y,m){
  if(any(is.na(c(y,m)))){
    NA
  } else if(m<10){
    glue("0{m}-{y}")
  } else {
    glue("{m}-{y}")
  }
}

prepr.tbl<-data.preprocessing(input.tbl)
save(prepr.tbl, file = "prepr.tbl.rda")


#PART 2###########################
input.tbl<- prepr.tbl
#' @keywords internal
#' @export exclud.sympt.comorb.tret
exclud.sympt.comorb.tret <- function(input.tbl){
  
  tot=nrow(input.tbl)
  tot_icu=nrow(filter(input.tbl,slider_icu_ever==TRUE))
  
  data<-dplyr::select(input.tbl, c(starts_with("symptoms_"),starts_with("comorbid_"),starts_with("treat_"))) %>%
    pivot_longer(c(starts_with("symptoms_"),starts_with("comorbid_"),starts_with("treat_")), 
                 names_to = "variable", 
                 values_to = "value")%>%
    mutate(count=1)%>%
    group_by(variable,value)%>%
    summarise(n = sum(count, na.rm=T))%>%
    mutate(prop=round(n/tot,digit=2))%>%
    filter(is.na(value))%>%
    filter(prop>=0.90)%>% 
    dplyr::select(variable)
  
  data2<-dplyr::select(input.tbl, c(starts_with("icu_treat"),slider_icu_ever)) %>%
    filter(slider_icu_ever==TRUE)%>%
    pivot_longer(c(starts_with("icu_treat")), 
                 names_to = "variable", 
                 values_to = "value")%>%
    mutate(count=1)%>%
    group_by(variable,value)%>%
    summarise(n = sum(count, na.rm=T))%>%
    mutate(prop=round(n/tot_icu,digit=2))%>%
    filter(is.na(value))%>%
    filter(prop>=0.90)%>%
    dplyr::select(variable)
  
  rmv<-unique(c(data$variable, data2$variable))
  
}

input.tbl<-input.tbl%>%
  dplyr::select(-c(all_of(rmv)))

save(input.tbl, file = "input.tbl_partner.rda")

