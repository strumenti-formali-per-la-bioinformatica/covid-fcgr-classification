library(readr)
library(dplyr)
library(ggplot2)
library(survival)
library(survminer)

clade_G <- read_csv("C:/Users/carme/OneDrive/Desktop/UNIVERSITA 2/BIOINF/clade_G.csv")
clade_G$Type<-as.factor(clade_G$Type)
clade_G$Clade<-as.factor(clade_G$Clade)
clade_G$`Pango lineage`<-as.factor(clade_G$`Pango lineage`)
clade_G$`Pango version`<-as.factor(clade_G$`Pango version`)
clade_G$`AA Substitutions`<-as.factor(clade_G$`AA Substitutions`)
clade_G$Variant <-as.factor(clade_G$Variant )
clade_G$Location<-as.factor(clade_G$Location)
clade_G$Gender<-as.factor(clade_G$Gender          )
clade_G$`Patient status`<-as.factor(clade_G$`Patient status`)
clade_G$`Passage details/history`<-as.factor(clade_G$`Passage details/history`)
clade_G$`Sampling strategy`<-as.factor(clade_G$`Sampling strategy`)
clade_G$Treatment<-as.factor(clade_G$Treatment)
clade_G$`Specimen source`<-as.factor(clade_G$`Specimen source`)


##ANALISI DEI MISSING VALUES 
sum(is.na(clade_G))
summary(clade_G)



## TEST DEL CHI QUADRO RISPETTO ALLA VARIANTE
table(clade_G$Variant)

variant_table <- table(clade_G$Variant)
sorted_variant_table <- sort(variant_table, decreasing = TRUE)
top_5_variants <- names(sorted_variant_table)[1:4]
filtered_clade_G <- clade_G %>% filter(Variant %in% top_5_variants)
print(top_5_variants)
filtered_clade_G$`Variant` <- factor(filtered_clade_G$`Variant`, levels = top_5_variants)


vcramer <- function (x, y = NULL) 
{
  if(!is.null(y)) {
    tab <- table(x, y)
  } else tab = as.matrix(x)
  
  n <- (min(nrow(tab), ncol(tab))-1) * margin.table(tab)
  chiq <- as.numeric(chisq.test(tab, correct = FALSE)$statistic)
  p <- chisq.test(tab, correct = FALSE)$p.value
  v = sqrt(chiq / n)
  res <- c("chi.sq" = chiq, "p" = p, "v di Cramer" = v)
  res
}



filtered_clade_G$Variant <- gsub("Former (VOC|VOI) ([A-Za-z]+).*", "\\2", filtered_clade_G$Variant)



#Con il test del chi quadrato di Pearson, si intende controllare se l'associazione
#fra due variabili, eventualmente evidenziata da una tabella di contingenza, sia 
#statisticamente significativa. Si tratta di un test di verifica di ipotesi che 
#attribuisce un valore di probabilità all'ipotesi nulla (cioè all'ipotesi di assenza di associazione).
pango_table <- table(filtered_clade_G$`Pango lineage`)
sorted_pango_table <- sort(pango_table, decreasing = TRUE)
top_5_pango <- names(sorted_pango_table)[1:3]
filtered_top_5_pango <- filtered_clade_G %>% filter(`Pango lineage` %in% top_5_pango)
print(top_5_pango)
print(filtered_top_5_pango)
filtered_top_5_pango$`Pango lineage` <- factor(filtered_top_5_pango$`Pango lineage`, levels = top_5_pango)


tab1<-table(filtered_top_5_pango$Variant, filtered_top_5_pango$`Pango lineage`)

round(prop.table(tab1, margin = 2)*100, 1)
res1 <- chisq.test(tab1)
res1




pango_table <- table(filtered_clade_G$`Pango version`)
sorted_pango_table <- sort(pango_table, decreasing = TRUE)
top_5_pango <- names(sorted_pango_table)[1:3]
filtered_top_5_pango <- filtered_clade_G %>% filter(`Pango version` %in% top_5_pango)
print(top_5_pango)
print(filtered_top_5_pango)
filtered_top_5_pango$`Pango version` <- factor(filtered_top_5_pango$`Pango version`, levels = top_5_pango)


tab1<-table(filtered_top_5_pango$Variant, filtered_top_5_pango$`Pango version`)

round(prop.table(tab1, margin = 2)*100, 1)
res1 <- chisq.test(tab1)
res1
tab <- table(filtered_top_5_pango$Variant,filtered_top_5_pango$`Pango version`)
tab1<-table(filtered_top_5_pango$`Pango version`)
prop.table(tab1)
tab
tab
prop.table(tab[1,])
tab
prop.table(tab[2,])
tab
prop.table(tab[3,])
tab
prop.table(tab[4,])





pango_table <- table(filtered_clade_G$`Gender`)
sorted_pango_table <- sort(pango_table, decreasing = TRUE)
top_5_pango <- names(sorted_pango_table)[2:3]
filtered_top_5_pango <- filtered_clade_G %>% filter(`Gender` %in% top_5_pango)
print(top_5_pango)
print(filtered_top_5_pango)
filtered_top_5_pango$`Gender` <- factor(filtered_top_5_pango$`Gender`, levels = top_5_pango)


tab1<-table(filtered_top_5_pango$Variant, filtered_top_5_pango$`Gender`)

round(prop.table(tab1, margin = 2)*100, 1)
res1 <- chisq.test(tab1)
res1
tab <- table(filtered_top_5_pango$Variant,filtered_top_5_pango$`Gender`)
tab1<-table(filtered_top_5_pango$`Gender`)
prop.table(tab1)
tab
tab
prop.table(tab[1,])
tab
prop.table(tab[2,])
tab
prop.table(tab[3,])
tab
prop.table(tab[4,])






pango_table <- table(filtered_clade_G$`Patient status`)
sorted_pango_table <- sort(pango_table, decreasing = TRUE)
top_5_pango <- names(sorted_pango_table)[2:3]
filtered_top_5_pango <- filtered_clade_G %>% filter(`Patient status` %in% top_5_pango)
print(top_5_pango)
print(filtered_top_5_pango)
filtered_top_5_pango$`Patient status` <- factor(filtered_top_5_pango$`Patient status`, levels = top_5_pango)


tab1<-table(filtered_top_5_pango$Variant, filtered_top_5_pango$`Patient status`)

round(prop.table(tab1, margin = 2)*100, 1)
res1 <- chisq.test(tab1)
res1
tab <- table(filtered_top_5_pango$Variant, filtered_top_5_pango$`Patient status`)
tab1<-table(filtered_top_5_pango$`Patient status`)
prop.table(tab1)
tab
tab
prop.table(tab[1,])
tab
prop.table(tab[2,])
tab
prop.table(tab[3,])
tab
prop.table(tab[4,])



filtered_clade_G$`Patient age`[filtered_clade_G$`Patient age` == "unknown"] <- NA
filtered_clade_G <- filtered_clade_G %>% filter(!is.na(`Patient age`))
filtered_clade_G$`Patient age`<- as.numeric(filtered_clade_G$`Patient age`)
summary(filtered_clade_G$`Patient age`)
filtered_clade_G$age_category <- cut(filtered_clade_G$`Patient age`,
                                    breaks = c(0, 10, 28,50, Inf),
                                    labels = c("0-10", "11-28", "29-50",  "50+"),
                                    right = FALSE)
table(filtered_clade_G$age_category)
age_variant_table <- table(filtered_clade_G$Variant, filtered_clade_G$age_category)
tab1<-table(filtered_clade_G$Variant, filtered_clade_G$age_category)
round(prop.table(tab1, margin = 2)*100, 1)
res1 <- chisq.test(tab1)
res1
age_variant_table <- table(filtered_clade_G$Variant, filtered_clade_G$age_category)
age_variant_table
prop.table(age_variant_table[1,])
age_variant_table
prop.table(age_variant_table[2,])
age_variant_table
prop.table(age_variant_table[3,])
age_variant_table
prop.table(age_variant_table[4,])











pango_table <- table(filtered_clade_G$`Specimen source`)
sorted_pango_table <- sort(pango_table, decreasing = TRUE)
top_5_pango <- names(sorted_pango_table)[1:3]
filtered_top_5_pango <- filtered_clade_G %>% filter(`Specimen source` %in% top_5_pango)
print(top_5_pango)
print(filtered_top_5_pango)
filtered_top_5_pango$`Specimen source` <- factor(filtered_top_5_pango$`Specimen source`, levels = top_5_pango)


tab1<-table(filtered_top_5_pango$Variant, filtered_top_5_pango$`Specimen source`)

round(prop.table(tab1, margin = 2)*100, 1)
res1 <- chisq.test(tab1)
res1
tab <- table(filtered_top_5_pango$Variant, filtered_top_5_pango$`Specimen source`)
tab1<-table(filtered_top_5_pango$`Specimen source`)
prop.table(tab1)
tab
tab
prop.table(tab[1,])
tab
prop.table(tab[2,])
tab
prop.table(tab[3,])
tab
prop.table(tab[4,])

