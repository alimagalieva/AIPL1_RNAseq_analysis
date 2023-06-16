
rm(list=ls(all=TRUE))

library(ggplot2)
library(tidyverse)
library(plyr) 
library(dplyr)
library(ggrepel)
library(stringr)  
library(data.table)

df <- read.csv('Clinical Trials Updated feb2023_K.csv', sep = ';') 
df <- df[-c(263:269),]

#indication_sum <- data.frame(table(df$Indication))
#colnames(indication_sum)[1] ="Indication"

indication_map <- read.csv('41434_2022_363_MOESM4_ESM.csv', sep = ';') 

df <- merge(df, indication_map, 
            by = 'Indication', 
            all.x = T)

area <- data.frame(table(df$Therapeutic.Area))

adm <- data.frame(table(df$Administration))

serotype <- data.frame(table(df$AAV.serotype))

promoter <- data.frame(table(df$Promoter))

safety <- data.frame(table(df$Safety.met))

efficacy <- data.frame(table(df$Primary.Efficacy.Endpoint.Met))












################################################################################
############################### generating table ###############################
################################################################################

rm(list=ls(all=TRUE))

library(ggplot2)
library(tidyverse)
library(plyr) 
library(dplyr)
library(ggrepel)
library(stringr)  
library(data.table)

df2023 <- read.csv('SearchResults_AAV_2023.csv', sep = ',') 
df2022 <- read.csv('41434_2022_363_MOESM3_ESM_2022.csv', sep = ';', 
                   stringsAsFactors=FALSE, fileEncoding="latin1")
df2023 <- df2023[,!(names(df2023) %in% c('Rank'))]
colnames(df2023)[1] ="Clinical.Trial.Identifier"

df <- merge(df2023, df2022, 
            by = 'Clinical.Trial.Identifier', 
            all = T)
#df <- df[,!(names(df) %in% c('Acronym', 'Sponsor.Collaborator', 'Other.IDs', 'Study.Documents', 'URL', 'Conditions', 'Status.y', 'Primary.Completion.Date.y'))]

df_filt <- filter(df, Study.Type %in% c('Interventional', NA))

## tried to delete automatically rows with ANCA-ass.vasc, but never did it
#drop <- df_filt[str_detect(df_filt, 'Associated Vasculitis'), ]
# or: 
# df_filt <- df_filt[!df_filt %like% 'Associated Vasculitis', ] 

#df_filt = df_filt[(!row.names(df_filt) %in% row.names(drop)), ] 

write.table(df_filt, file="data2023_filt.txt", sep="\t", quote=F) 

## added info manually, saved to new file 


################################################################################
#################################### Fig. 5 ####################################
################################################################################

df <- data.frame(condition = c('musculoskeletal', 'neurology', 'cardiovascular', 
                               'haematology', 'metabolic', 'ophthalmology', 
                               'other', 'virology', 'oncology', 'autoimmune'), 
                 number = c(10, 17, 4, 14, 20, 26, 1, 2, 2, 3))

df <- df %>% 
  arrange(desc(condition)) %>%
  mutate(prop = number / sum(df$number)*100) 

#jpeg("5_plot.jpg", res=100)

df %>% 
  ggplot(aes(x = "", y = prop, fill = reorder(condition, prop))) +
   geom_bar(stat="identity", width=2, color="white") +
   coord_polar(theta = "y") + 
   geom_text(aes(label = paste0(number, '%'), x=2.3), 
             position = position_stack(vjust = 0.5), 
             color = "black", size=4) +
   labs(x = NULL, y = NULL, fill = NULL, title = "") +
   guides(fill = guide_legend(reverse = TRUE)) +
   theme_void()  # remove background, grid, numeric labels

#dev.off()

'''
df <- df %>% 
  arrange(desc(condition)) %>%
  mutate(prop = number / sum(df$number)*100) %>%
  mutate(ypos = cumsum(number)- 0.5*prop)
  
ggplot(df[1:10,], aes(x = "", y = prop, fill = condition)) +
  geom_bar(stat="identity", width=2, color="white") +
  coord_polar("y", start = 0) +
  geom_text(aes(y = ypos, x=2.3, label = paste0(number, '%'), color = "black", 
  size=4)) +
  labs(x = NULL, y = NULL, fill = NULL, title = "") +
  theme_void()  
'''

################################################################################
################################# Fig. 1 donut #################################
################################################################################

rm(list=ls(all=TRUE))

df <- data.frame(vector = c('AAV', 'AdV', 'HSV', 'others', 'lentivirus'), 
                 number = c(137, 83, 46, 12, 6))

df <- df %>% mutate(prop = round(number/sum(df$number)*100, 1)) 

# Compute the cumulative percentages (top of each rectangle)
df$ymax <- cumsum(df$prop)
# Compute the bottom of each rectangle
df$ymin <- c(0, head(df$ymax, n=-1))
# Compute label position
df$labelPosition <- (df$ymax + df$ymin) / 2

df %>% ggplot(aes(x = 5, y = prop, fill = reorder(vector, prop))) +
  geom_col(width=3.5) +
  xlim(0, 10.5) +
  scale_fill_viridis_d() +
  coord_polar(theta = "y") +
  guides(fill = guide_legend(reverse = TRUE)) +
  labs(x = NULL, y = NULL, fill = NULL, title = "") +
  geom_text_repel(aes(label = paste0(prop, '%\n', vector)), 
                  color="black", size = 4, x = 10.5, y = df$labelPosition) +
  geom_text(aes(label = ' Vector \n type '),
            color = 'black', size = 8, x = 0, y = 0) +
  theme_void()
# and then move text in Paint 


################################################################################
#################################  Fig. 1 Hist #################################

ggplot(df, aes(x=vector, y=number, fill = reorder(vector, prop))) +
  geom_bar(stat="identity") +
  geom_text(aes(label=number), vjust=-0.5, color="black",
            position = position_dodge(0.9), size=4.5) +
  ggtitle("N = 284 clinical trials") +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("vector type") + 
  ylab("# of active clinical trials") +
  guides(fill = guide_legend(reverse = TRUE)) +
  scale_fill_viridis_d(option = "D", aesthetics = "fill", begin = 0, end = 1) +
  theme_minimal()
  

################################################################################
#################################### Fig. 7 ####################################
################################################################################

rm(list=ls(all=TRUE))

df <- data.frame(administration = c('others', 'intracoronary', 'intraarticular', 
                                    'intrathecal', 'intravitreal', 'intramuscular', 
                                    'subretinal', 'intracraneal', 'intravenous'), 
                 prop = c(2, 3, 4, 4, 9, 9, 17, 17, 35))

df %>% 
  ggplot(aes(x = "", y = prop, fill = reorder(administration, prop))) +
  geom_bar(stat="identity", width=2, color="white") +
  coord_polar(theta = "y") + 
  geom_text(aes(label = paste0(prop, '%'), x=2.3), 
            position = position_stack(vjust = 0.5), 
            color = "black", size=5) +
  labs(x = NULL, y = NULL, fill = NULL, title = "") +
  guides(fill = guide_legend(reverse = TRUE)) +
  theme_void()  # remove background, grid, numeric labels





