# Load all necessary packages  -----------------------------------------------------

pacman::p_load(easyPubMed,
               tidyverse,
               broom,
               countrycode,
               WDI,
               rio,
               data.table,
               ggrepel,
               psych,
               slopegraph,
               PerformanceAnalytics,
               maptools,
               ggthemr,
               knitr)
# End of command chunk ----------------------------------------------------------------


# Load world map & Africa map data --------------------------------------------------              
# To get Africa countries names and major alternative names

data(wrld_simpl)
afr=wrld_simpl[wrld_simpl$REGION==2,]
rm(wrld_simpl)

country.name <- data.frame(country = afr$NAME) %>%
  mutate(country_names = recode(country,
                                "Benin" = "Benin NOT Nigeria",
                                "Cape Verde" = "Cape Verde OR Cabo Verde",
                                "Cote d'Ivoire" = "Cote d'Ivoire OR 'ivory coast'",
                                "Niger" = "Niger NOT nigeria",
                                "Swaziland" = "Swaziland OR Swasiland OR Eswatini"),
         iso3c = countrycode(country, 'country.name', 'iso3c'))

country <- as.character(country.name$country_name)
# End of command chunk --------------------------------------------------------------------


#  Main variables declaration  -------------------------------------------------------------  
topic <- "primary health care [MeSH] OR primary health care"          # Topic of Interest
year <- 1999:2019                       # Study period
indicators <- c("SP.POP.TOTL",          # Country's population     
                "NY.GDP.MKTP.CD")       # Country's GDP  
# End of command chunk --------------------------------------------------------------------



# Get world bank indicator data for Pop & GDP  ---------------------------------------------
wdi_data1 <- WDI(indicator = indicators, start = 1999, end = 2019, extra = TRUE) %>%
  filter(region  != "Aggregates") %>% 
  mutate(continent = countrycode(iso3c, 'iso3c', 'continent')) %>% 
  filter(continent == "Africa") %>% select(-c(continent)) %>%
  rename(pop =  SP.POP.TOTL, 
         gdp = NY.GDP.MKTP.CD) %>% arrange(country, year) %>% select(-country) 
# End of command chunk --------------------------------------------------------------------




# Define Function to collect data for Worldwide publication on the topic, yearly ------------
yearFunc <- function(year) { 
  query <- paste0(topic, " AND ", year, "[PDAT]")
  fetch <- get_pubmed_ids(query)
  
  tibble(year = year,
         publ_all = as.numeric(fetch$Count)) 
} 

df_world <- map_df(year, yearFunc)
# End of command chunk --------------------------------------------------------------------



# Main function for topic yearly data collection for each country ---------------------------
biblioFunc <- function(country) {
  
  map_df(year, function(year) { 
    query <- paste0(topic, " AND ", country, " AND ", year, "[PDAT]")
    fetch <- get_pubmed_ids(query)
    
    tibble(country_names = country,
           year = year,
           publ = as.numeric(fetch$Count)) } )
}

df_afr <-  map_df(country, biblioFunc)
# End of command chunk --------------------------------------------------------------------


# Further data cleaning, merge with WDI data ------------------------------------------------

df_afr <- df_afr %>%
  left_join(., country.name) %>% select(country, everything(), -country_names) %>%
  mutate(iso3c = countrycode(country, 'country.name', 'iso3c'),
         continent = countrycode(iso3c, 'iso3c', 'continent'),
         georegion = countrycode(iso3c, 'iso3c', 'region'),
         yearm = year - 1998,
         period = cut(year, breaks = 4, include.lowest = TRUE, 
                      labels = c("1999-2003", "2004-2008", "2009-2013", "2014-2019"))
         ) %>% 
  left_join(., wdi_data, by = c("iso3c", "year"))  %>%
  fill(everything())
# End of command chunk --------------------------------------------------------------------



# Generate main data per country ------------------------------------------------------------
df_main <- df_afr %>%
  group_by(country) %>%
  summarise(
    georegion = first(georegion),                                    
    iso2c = first(iso2c),
    iso3c = first(iso3c),
    total_publ = sum(publ),                                         #Sum publication for study period
    start_publ = first(publ),                                       #Start No of publication
    end_publ = last(publ),                                          #End No of publication
    abs_incr = (end_publ / start_publ)*100,                         #Absolute increase
    rel_incr = ((end_publ - start_publ)/ start_publ)*100,           #Relative increase
    pop = last(pop),
    gdp = last(gdp),
    income = first(income),
    adj_pop = (total_publ / pop) * 1000000,                         #Adjusted for population
    adj_gdp = (total_publ / gdp) * 1000000000,                      #Adjusted for GDP
    pub_cat = cut(total_publ,                                       #Informative quantile for No of Publication
                       breaks = c(0, 10, 50, 100, 500, 20000),
                       labels = c("1 to 10", "11 to 50", "51 to 100", 
                                  "101 to 500", "501 to 3500")))  %>% arrange(total_publ)
# End of command chunk --------------------------------------------------------------------



# Collate the list of top countries, absolute number, adj for POP, adj for GDP
top10     <- df_main %>% arrange(desc(total_publ)) %>% slice(1:10) %>% .$iso3c
top10_pop <- df_main %>% arrange(desc(adj_pop)) %>% slice(1:10) %>% .$iso3c
top10_gdp <- df_main %>% arrange(desc(adj_gdp)) %>% slice(1:10) %>% .$iso3c
# End of command chunk --------------------------------------------------------------------



# Calculate AAPC - average annual percentage change
# From poisson regression
# Nested data ------------------------------------------------------------------------------

df_nest <- df_afr %>% filter(iso3c %in% top10) %>%
  group_by(country, georegion) %>%
  nest()
# End of command chunk --------------------------------------------------------------------


# Fit models -------------------------------------------------------------------------------
country_model <- function(df) {
  glm(publ ~ yearm , family = poisson, data = df)
}


models <- df_nest %>%
  mutate(
    model  = data %>% map(country_model),
    tidy    = model %>% map(broom::tidy, exponentiate = TRUE, conf.int = TRUE))
  
df_aapc <- unnest(models, tidy) %>% filter(term == "yearm") %>%
  mutate(pvalue =sprintf("%4.3f", p.value),
         aapc = 100*(estimate -1),
         aapc_effect = paste0(sprintf("%4.1f", 100*(estimate-1)), " (", 
                              sprintf("%4.1f", 100*(conf.low-1)), " to ", 
                              sprintf("%4.1f", 100*(conf.high-1)), ")")) %>%
  select(country, aapc, aapc_effect, pvalue) 

# End of command chunk --------------------------------------------------------------------



# Export Table 1 --------------------------------------------------------------------------
df_table1 <- df_main %>% 
  filter(iso3c %in% top10)  %>% 
  left_join(., df_aapc) %>%
  select(country, total_publ, start_publ, end_publ, abs_incr, rel_incr, aapc_effect, pvalue)

kable(df_table1)
write.table(df_table1, file = "Ftable1.txt", sep = ",", quote = FALSE, row.names = F)
# End of command chunk --------------------------------------------------------------------

data(wrld_simpl)
afr=wrld_simpl[wrld_simpl$REGION==2, 142,]
rm(wrld_simpl)

test2 <- as_tibble(wrld_simpl@data)
test3 <- afr@data

test4 <- test2 %>% left_join(., wdi_data1, by = c('ISO3' = 'iso3c')) %>%

wdi_data1 <- WDI(indicator = indicators, start = 1998, end = 2017, extra = TRUE) %>%
  filter(region  != "Aggregates") %>% 
  mutate(continent = countrycode(iso3c, 'iso3c', 'continent')) %>% 
  filter(continent == "Africa") %>% select(-c(continent)) %>%

# Generate Map of Africa  - Figure 1 ------------------------------------------------------
afr@data <- afr@data %>% left_join(., df_main, by = c('ISO3' = 'iso3c'))
df_map  <- tidy(afr) 
df_map <- df_map %>% left_join(. , afr@data, by=c("id"="ISO3"))
cnames <- aggregate(cbind(long, lat) ~ NAME, data=df_map, FUN=function(x) mean(range(x))) 

ggplot() +
  geom_polygon(data = df_map, 
               aes(x = long, y = lat, group = group, fill = pub_cat.y),
               color = "black", size = 0.2) +
  scale_fill_manual(values = c("gray", "yellow", "orange", "pink", "blue")) + 
  coord_equal() +
  labs(x = NULL, y = NULL, fill = "No of Publications") +
  geom_text_repel(aes(x = long, y = lat, label =NAME), data = cnames, size = 2.5) + 
  theme_void() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        text = element_text(size=15),
        legend.position = c(0.9, 0.2)) 

# End of command chunk --------------------------------------------------------------------


# Export the figure 
dev.copy2pdf(file="Fig1_map.pdf")
dev.off()
# End of command chunk --------------------------------------------------------------------




# GENerate Trend Plot overall publications & percent of worldshare -----------------------
df_trend <- df_afr %>%
  group_by(year) %>%
  summarise(publ_afr = sum(publ),
            period = first(period)) %>%
  left_join(., df_world) %>%
  mutate(perc = (publ_afr / publ_all)*100,
         yearm = year - 1968,
         yearm2 = yearm^2)

scaleFactor = (max(df_trend$publ_afr) / max(df_trend$perc))
# End of command chunk --------------------------------------------------------------------

# Generate the plot
ggthemr('fresh')
df_trend %>%
  ggplot(., aes(x = year)) +
  geom_bar(aes(y = publ_afr), stat = "identity") +
  geom_point(aes(y = publ_afr), size = 6, shape = 21, fill = "white") +
  geom_line(aes(y = perc*scaleFactor), color = "red", size = 3)  +
  geom_point(aes(y = perc*scaleFactor), size = 6, shape = 21, fill = "white") +
  geom_text(aes(y = perc*scaleFactor, label=sprintf("%4.1f", perc), hjust=0, vjust=3.5)) +
  scale_y_continuous(name = "Total publication", 
                     sec.axis = sec_axis(~./scaleFactor, name = "% World share")) +
  scale_x_continuous("Publication year", breaks = 1999:2019) +
  theme(axis.line.y.right = element_line(color = "red"), 
        axis.ticks.y.right = element_line(color = "red"),
        axis.text.y.right = element_text(color = "red"), 
        axis.title.y.right = element_text(color = "red"),
        text = element_text(face="bold", size=20),
        axis.text.x = element_text(face="bold", color="black", size=18, angle=45),
        axis.text.y = element_text(face="bold", color="black", size=18))
# End of command chunk --------------------------------------------------------------------


# Export the figure 
dev.copy2pdf(file="Fig2_trend.pdf")
dev.off()
# End of command chunk --------------------------------------------------------------------




# Generate Trends fig for different period -------------------------------------------------
df_period <- df_afr %>% 
  group_by(country, period) %>%
  summarise(publ = sum(publ)) %>% 
  dcast(., country  ~ period, value.var = "publ") %>%
  mutate_if(is.numeric, funs(rank(desc(.), ties.method = 'first'))) %>%
  arrange((`1999-2003`)) %>%
  as.data.frame(.) %>%
  column_to_rownames(var = 'country') 

cols <- `[<-`(rep("black", 57), 5, "red")
slopegraph(df_period, xlim = c(-1, 7), ylim = c(57,0), offset.x = 0.06,
           col.lines = cols, col.lab = cols, 
           main = 'Relative Rank of Primary Health Care Research in Africa, 1999 - 2019')
# End of command chunk --------------------------------------------------------------------


# export figure
dev.copy2pdf(file="Fig3_trend.pdf")
dev.off()
# End of command chunk --------------------------------------------------------------------






# Correlations with GDP and Population ----------------------------------------------------
# Correlation analysis
df_corr <- df_main %>% 
  select(total_publ, gdp,pop) %>%
  mutate_if(is.numeric, funs(log)) 

x <- as_vector(df_corr[,1])
y.gdp <-  as_vector(df_corr[,2])
y.pop <-  as_vector(df_corr[,3])


cor.test(x, y.gdp, method=c("pearson", "kendall", "spearman"))
cor.test(x, y.pop, method=c("pearson", "kendall", "spearman"))
# End of command chunk --------------------------------------------------------------------


# Quick plot
chart.Correlation(df_corr, pch=21)
# End of command chunk --------------------------------------------------------------------



# GENerate figure for Top 10 adjusted POP ----------------------------------------------------
ggthemr_reset()
df_afr %>% 
  mutate(adj = (publ / pop)*1000000) %>%
  filter(iso3c %in% top10_pop) %>%
  ggplot(., aes(x = fct_reorder(country, adj, .fun = median, .desc =TRUE), y = adj)) + 
  geom_boxplot(aes(fill = fct_reorder(country, adj, .fun = median, .desc =TRUE))) + 
  geom_jitter(position=position_jitter(0.2)) +
  theme_bw(base_size = 14) +
  xlab("Countries") +
  ylab("Number of publications per 1,000, 000 population") +
  scale_fill_discrete(guide = guide_legend(title = "Country")) +
  theme(axis.text.x  = element_text(angle=10, vjust=0.5, size=16)) +
  guides(fill = FALSE)
# End of command chunk --------------------------------------------------------------------


# Export the figure
dev.copy2pdf(file="Fig4_topPOP.pdf")
dev.off()
# End of command chunk --------------------------------------------------------------------




#  GENerate Scatter plot for the correlation between POP & Publ  ----------------------------
df_main %>%
  ggplot(aes(x= log(total_publ), y = log(pop))) +
  geom_point(aes(size = gdp, color = georegion)) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_size(range = c(0, 10)) +
  theme(legend.position='bottom') +
  labs(x = "log(Country's population)",
       y ="log(Total publications)",
       size = "Sample size",
       color = "Region")  +
  theme(text = element_text(size=14))  +
  theme(axis.title=element_text(size=20,face="bold"))  +
  geom_text_repel(aes(label = iso3c),
                  box.padding = unit(0.2, "lines"))  +
  guides(size = FALSE)
# End of command chunk --------------------------------------------------------------------



# Export the figure
dev.copy2pdf(file="Fig5_corr.pdf")
dev.off()
# End of command chunk --------------------------------------------------------------------



# GENerate figure for Top 10 adjusted GDP ---------------------------------------------------
ggthemr_reset()
df_afr %>% 
  mutate(adj = (publ / gdp)*1000000) %>%
  filter(iso3c %in% top10_pop) %>%
  ggplot(., aes(x = fct_reorder(country, adj, .fun = median, .desc =TRUE), y = adj)) + 
  geom_boxplot(aes(fill = fct_reorder(country, adj, .fun = median, .desc =TRUE))) + 
  geom_jitter(position=position_jitter(0.2)) +
  theme_bw(base_size = 14) +
  xlab("Countries") +
  ylab("Number of publications per 1,000, 000 population") +
  scale_fill_discrete(guide = guide_legend(title = "Country")) +
  theme(axis.text.x  = element_text(angle=10, vjust=0.5, size=16)) +
  guides(fill = FALSE)
# End of command chunk --------------------------------------------------------------------


# export the figure
dev.copy2pdf(file="Fig6_topGDP.pdf")
dev.off()
# End of command chunk --------------------------------------------------------------------


# GENerate Scatter plot for the correlation between GDP & Publ -----------------------------
df_main %>%
  ggplot(aes(x= log(total_publ), y = log(gdp))) +
  geom_point(aes(size = gdp, color = georegion)) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_size(range = c(0, 10)) +
  theme(legend.position='bottom') +
  labs(x = "log(Gross Domestic Product)",
       y ="log(Total publications)",
       size = "Sample size",
       color = "Region")  +
  theme(text = element_text(size=14))  +
  theme(axis.title=element_text(size=20,face="bold"))  +
  geom_text_repel(aes(label = iso3c),
                  box.padding = unit(0.2, "lines"))  +
  guides(size = FALSE)
# End of command chunk --------------------------------------------------------------------


# export the figure
dev.copy2pdf(file="Fig7_gdp.pdf")
dev.off()

# End of command chunk --------------------------------------------------------------------










