################################################################################################################################################################
######Identification of top candidate windows and testing window-by_window signatures of convergence (null-W test) with considering recombination effect:#######
################################################################################################################################################################


#################################################################Packages#######################################################################################
library(tidyverse)
library(data.table)
library(foreach)
library(doParallel)
##################################################################Arguments#########################################################################################
args = commandArgs(trailingOnly = TRUE)

# argument1: association result (GWAS,BYPASS, or correlation) in species 1, includes for columns (with out header) as follow: 1. species  2. Variable 3.SNP_ID (example Ha412HOChr01:180057) 4.window ID (example Ha412HOChr01:180001-185000) 5.asociation value (P,BF or r)
# argument2: association result (GWAS,BYPASS, or correlation) in species 2, includes for columns (with out header) as follow: 1. species  2. Variable 3.SNP_ID (example Ha412HOChr01:180057) 4.window ID (example Ha412HOChr01:180001-185000) 5.asociation value (P,BF or r)
# argument3: association value , P for P value, R for correlation, BF for Bayes Factor

#example: Rscript LDCLUST.R annuus_GWAS_result.txt argophyllus_GWAS_result.txt P

##################################################################Files#########################################################################################

recomb <- "/data/users/mjahani/temp/5kbwindow_recombination" #recombination ratio

################################################################Theresholds######################################################################################

association_cut <- 0.01 #The threshold for the significant association. Example: for P-value, the 0.01 cut off means the bottom one percent quantile of the association tests (smaller P-value = stronger association), while for correlation, the threshold should be 0.99 (as larger correlation means stronger association)
dbinom_cut <- 0.9999 #The threshold for identifying top candidate windows (0.9999 quantile of the binomial expectation)
q <- 0.05 # the threshold for q-value in null-W results (Empirical p-values of null-W were converted to q-values for estimation of False Discovery Rate)
nthreads = 48 #Number of threads for parallel computation 

################################################################Output######################################################################################

#$argument1_topcandidate #the file indicates if each window is a top candidate for the species in argument 1, and contains information on number of SNPs (snp_count) and number of significantly associated snps(outlier_count) for each window
#$argument2_topcandidate #the file indicates if each window is a top candidate for the species in argument 2, and contains information on number of SNPs (snp_count) and number of significantly associated snps(outlier_count) for each window
#$species1_$species2_$variablename_nullw_result #null-W test results for species 2 topcandidates in species 1 genome
#$species1_$species2_$variablename_nullw_result_FDR_0.05 #null-W test results for species 2 topcandidates in species 1 genome which passed the FDR test
#$species2_$species_1$variablename_nullw_result #null-W test results for species 1 topcandidates in species 2 genome
#$species2_$species1_$variablename_nullw_result_FDR_0.05 #null-W test results for species 1 topcandidates in species 2 genome which passed the FDR test

#################################################Identification of top candidate windows######################################################################################

registerDoParallel(cores=nthreads)
data <- foreach(i=1:2, .combine='rbind', .errorhandling='stop') %dopar% {
fread(args[i],
      showProgress = F,
      header = F) %>% 
  select(species = V1,
         variable = V2,
         window = V4,
         value = V5) %>%
  mutate(outlier = ifelse(value <= quantile(as.numeric(value) ,
                                            association_cut)  ,
                  "outlier" ,
                  "non_outlier")) %>%
  group_by(species,
           variable,
           window,
           outlier) %>% 
  summarise(snp_count = n()) %>%
  ungroup() %>%
  spread(outlier,
         snp_count) %>%
  mutate(outlier = replace_na(outlier, 0) ,
         non_outlier = replace_na(non_outlier,
                                  0)) %>%
  mutate(snp_count = as.numeric(outlier) + as.numeric(non_outlier)) %>%
  select(-non_outlier) %>%
  rename(outlier_count = outlier) %>%
  left_join(. , 
            rename(fread(recomb,sep="\t", header = F) , 
                   window = V1,
                   recomb_rate = V2)) %>%
  filter(!is.na(recomb_rate)) %>%
  mutate(recomb_quantile = cut(recomb_rate, 
                      breaks = quantile(recomb_rate ,
                                        probs=seq(0 ,
                                                  1 , 
                                                  by=0.2) ,
                                        na.rm=TRUE) , 
                      labels = c("0-20" ,
                                 "20-40" ,
                                 "40-60" ,
                                 "60-80" ,
                                 "80-100"), 
                      include.lowest=TRUE)) %>%
  select(-recomb_rate) -> SNP_count


snp_count2 <- left_join(SNP_count ,
                       SNP_count %>% 
                         group_by(species ,
                                  variable ,
                                  recomb_quantile) %>% 
                         summarise(totsnp = sum(snp_count) ,
                                   totout = sum(outlier_count)) %>% 
                         mutate(expected = totout/totsnp) %>% 
                         select(-totsnp,
                                -totout)) %>% 
  mutate(Qbinom = qbinom (dbinom_cut ,
                          snp_count ,
                          expected)) %>%
  mutate(top_candidate = ifelse(outlier_count > Qbinom,
                                "top",
                                "not")) %>%
  select(species,
         variable,
         window,
         outlier_count,
         snp_count,
         Qbinom,
         top_candidate) %>%
  fwrite(paste0(args[i],"_topcandidate"),
         sep = "\t",
         col.names = T)

rm(SNP_count)

fread(paste0(args[i],"_topcandidate"),
      showProgress = F,
      header = T) 
}

#######################################################################data preperation for null-w test################################################################

data %>%
  distinct(species) %>%
  pull() -> species_pair

rbind(data %>%
        group_by(window) %>%
        tally() %>%
        filter(n == 2) %>%
        select(window) %>%
        mutate(species = species_pair[1]) %>%
        mutate(species2 = species_pair[2]),
      data %>%
        group_by(window) %>%
        tally() %>%
        filter(n == 2) %>%
        select(window) %>%
        mutate(species = species_pair[2]) %>%
        mutate(species2 = species_pair[1])) -> orth_window_species 


data %>% 
  select(species2 = species ,
         variable ,
         window ,
         species2_top_candidate = top_candidate) -> top_compare

ortho_species1.2_top <-foreach(i=1:2, .combine='rbind', .errorhandling='stop') %dopar% {
  fread(args[i],
        showProgress = F,
        header = F) %>% 
    select(species = V1 ,
           variable = V2 ,
           snp_ID = V3 ,
           window = V4 ,
           value = V5) %>%
    left_join(.,
              data) %>%
    filter(!is.na(snp_count)) %>%
    select(-outlier_count ,
           -Qbinom) %>%
    left_join(.,
              orth_window_species) %>%
    filter(!is.na(species2)) %>%
    left_join(.,
              top_compare) %>%
    filter(!is.na(species2_top_candidate)) 
  }

rm(data,orth_window_species,top_compare)

for (i in c(1,2)) {

ortho_species1.2_top %>%
  filter(species == species_pair[i]) %>%
  select(window) %>%
  distinct(window) %>%
  left_join(.,
  rename(
    fread(recomb ,
          sep="\t" ,
          header = F) , 
         window = V1,
         recomb_rate = V2)) %>%
filter(!is.na(recomb_rate)) %>%
  mutate(recomb_quantile = cut(recomb_rate, 
                               breaks = quantile(recomb_rate ,
                                                 probs=seq(0 ,
                                                           1 , 
                                                           by=0.2) ,
                                                 na.rm=TRUE) , 
                               labels = c("0-20" ,
                                          "20-40" ,
                                          "40-60" ,
                                          "60-80" ,
                                          "80-100"), 
                               include.lowest=TRUE)) -> window_recom_ortho
ortho_species1.2_top %>%
  filter(species == species_pair[i]) %>%
select(variable,value,snp_ID,window,snp_count,top_candidate,species2_top_candidate) %>% 
  left_join(.,window_recom_ortho) %>%
  filter(!is.na(recomb_rate)) %>%
  select(-recomb_rate) -> two_variable_orth

rm(window_recom_ortho)

two_variable_orth %>%
  distinct(recomb_quantile) %>%
  pull(recomb_quantile) -> quantile_list
#################################################################################null-w test################################################################

Final_result <- NULL

for (z in 1:length(quantile_list)) {

  top_list <- two_variable_orth %>% 
    filter(recomb_quantile == quantile_list[z]) %>%
    filter(top_candidate=="top") %>% 
    distinct(window) %>%
    pull(window)

  background_10_ID <- two_variable_orth %>%
    filter(recomb_quantile == quantile_list[z]) %>%
    filter(!window %in% top_list) %>%
    distinct(snp_ID) %>%
    sample_n(10000) %>%
    pull(snp_ID)
  
  background_10KB <- two_variable_orth %>%
    filter(recomb_quantile == quantile_list[z]) %>%
    filter(snp_ID %in% background_10_ID) %>%
    pull(value) 
  
  window_list <- two_variable_orth %>%
    filter(recomb_quantile == quantile_list[z]) %>%
    filter(!window %in% top_list) %>%
    filter(!is.na(value)) %>%
    filter(snp_count > 5) %>%
    distinct(window) %>%
    pull(window)

  Values<-two_variable_orth %>%
    filter(recomb_quantile == quantile_list[z]) %>%
    select(window,value)

  window_list1 <- two_variable_orth %>%
    filter(recomb_quantile == quantile_list[z]) %>%
    filter(species2_top_candidate=="top") %>%
    filter(snp_count > 3) %>%
    distinct(window) %>%
    pull(window)

  N_window <- two_variable_orth %>%
    filter(recomb_quantile == quantile_list[z]) %>%
    filter(species2_top_candidate == "top") %>%
    distinct(window) %>%
    nrow()
  #data frame of windows and log values for first vriable   
  Values1 <- two_variable_orth %>%
    filter(recomb_quantile == quantile_list[z]) %>%
    select(window,value)
  
  if (args[3] == "P") {
    
    background_10KB %>%
      log10(.)*-1 -> background_10KB
    
    Values %>%
      mutate(Value=log10(value)*-1) -> Values
    
    Values1 %>%
      mutate(value=log10(value)*-1) -> Values1 
    
  } else if (args[3] == "R") { 
    
    background_10KB %>%
      mutate(Value=Value^2) -> background_10KB
    
    Values %>%
      mutate(Value=Value^2) -> Values
    
    Values1 %>%
      mutate(Value=Value^2) -> Values1 
    
  } else if (args[3] == "PF") { 
    
  } else { 
    print("Wrong value in the third argument")
    }
  
  if (length(window_list1) > 0) { 
  
    nulldist<-foreach(k = 1:length(window_list), .combine='rbind', .errorhandling='stop') %dopar% {#null distribution
      
      data.frame(N= sum(is.na(pull(filter(Values,window==window_list[k])))==F),#how many snp for each not top window with  atleast  5 snps
                 W=as.vector(wilcox.test(background_10KB, pull(filter(Values,window==window_list[k])))$statistic)) %>% #wilcox test, x=(10000 p_values in first variable that are not in top windows for first species) y=(P_values of first non_top window in first variable)
        filter(!is.na(W)) %>%                                                                                           
        mutate(Zscore= (2*W-10000*N) / sqrt(10000*N*(10000+N+1)/3)) %>%
        select(Zscore)
    }#null distribution
    
    rm(Values,window_list)
    
    Result <- foreach(l=1:length(window_list1), .combine='rbind', .errorhandling='stop') %dopar% {
      
      data.frame(species1=species_pair[i],
                 species2=species_pair[species_pair != species_pair[i]],
                 Variable=distinct(two_variable_orth,variable),
                 recomb_quantile=as.character(quantile_list[z]),
                 Window=window_list1[l],
                 W=as.vector(wilcox.test (pull(filter(Values1,window==window_list1[l])),background_10KB)$statistic),#wilcox test, x=(values of first not top window for the second species), y=(10000 p_values in first variable that are not in top windows for first species)
                 P=wilcox.test (pull(filter(Values1,window==window_list1[l])),background_10KB)$p.value,
                 Mean_test=mean(pull(filter(Values1,window==window_list1[l]))),
                 N_samples=length(pull(filter(Values1,window==window_list1[l]))[!is.na(pull(filter(Values1,window==window_list1[l])))]),
                 Mean_BG=mean(background_10KB),
                 N_BG=length(background_10KB[!is.na(background_10KB)]),
                 Count=N_window)
    }
    
    temporary_result <- Result %>% mutate(Pred_P=(2* W - N_BG * N_samples)/sqrt(N_BG * N_samples*(N_BG+N_samples+1)/3))
    
    rm(Values1,window_list1,background_10KB)
    
    Emp_result <- foreach(m=1:nrow(temporary_result), .combine='rbind', .errorhandling='stop') %dopar% {
      data.frame(Emp_p = 1-(sum(as.numeric(temporary_result[m,13]) > pull(nulldist,Zscore))/nrow(nulldist)),Window=temporary_result[m,5])
      }
    
    left_join(temporary_result
              ,Emp_result) %>%
      filter(!is.na(Emp_p)) %>%
      rbind(. ,
            Final_result) -> Final_result

    rm(nulldist,Result,temporary_result,Emp_result,N_window)} else {
    data.frame(variable = distinct(two_variable_orth,variable) ,
               recomb_quantile = quantile_list[z] ,
               species1 = species_pair[i] ,
               species2 = species_pair[species_pair != species_pair[i]]) %>%
      fwrite(file = paste(species_pair[1] ,
                          species_pair[2] ,
                          distinct(two_variable_orth,variable) ,
                          "nullw_result_No_outlier" ,
                          sep = "_") ,
             append = T ,
             sep = "\t")
    rm(background_10KB,window_list,Values,window_list1,N_window,Values1)}
  
}#quantile loop 

#################################################################################q-value calcylation and FDR test################################################################

L = length(Final_result[,14])
Final_result_sorted <- Final_result[order(Final_result[,14]),]
w = Final_result_sorted[Final_result_sorted[,14] < q * ((1:L) / L),]
fwrite(Final_result ,
       file = paste(species_pair[i] ,
                    species_pair[species_pair != species_pair[i]] ,
                    distinct(two_variable_orth,variable) ,
                    "nullw_result" ,
                    sep = "_") ,
       append = T ,
       sep = "\t")

fwrite(w ,
       file = paste(species_pair[i] ,
                    species_pair[species_pair != species_pair[i]] ,
                    distinct(two_variable_orth,variable) ,
                    "nullw_result_FDR" ,
                    q ,
                    sep = "_") ,
       append = T ,
       sep = "\t")

rm(Final_result ,
   L ,
   Final_result_sorted ,
   w)
}#species loop

