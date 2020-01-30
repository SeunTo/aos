pacman::p_load(easyPubMed,
               tidyverse,
               broom,
               countrycode,
               WDI,
               rio,
               data.table,
               ggrepel,
               psych,
               stringr,
               slopegraph,
               PerformanceAnalytics,
               maptools,
               ggthemr,
               knitr)


map <- Dat1 %>%
  mutate(iso3c = countrycode(Country, 'country.name', 'iso3c'),
         Cat_cat = cut(NumCat,                                       #Informative quantile for No of Publication
                       breaks = c(0, 50000, 1000000, 10000000, 100000000, 500000000),
                       labels = c("1 to 50000", "51000 to 1000000", "1000001 to 10000000", 
                                  "10000001 to 100000000", "101000000 to 500000000")))
data(wrld_simpl)
map1 <- as_tibble(wrld_simpl@data) %>% inner_join(., map, by = c('ISO3' = 'iso3c'))
df_map  <- tidy(wrld_simpl)
df_maplmic <- df_map %>% inner_join(. , map1, by=c("id"="ISO3"))
cnames <- aggregate(cbind(long, lat) ~ id, data=df_maplmic, FUN=function(x) mean(range(x))) 

ggplot() +
  geom_polygon(data = df_maplmic, 
               aes(x = long, y = lat, group = group, fill = Cat_cat),
               color = "black", size = 0.2) +
  scale_fill_manual(values = c("gray", "yellow", "orange", "pink", "red")) + 
  coord_equal() +
  labs(x = NULL, y = NULL, fill = "NumCat") +
  geom_text_repel(aes(x = long, y = lat, label =id), data = cnames, size = 2.5) + 
  theme_void() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        text = element_text(size=15),
        legend.position = c(0.85, 0.1))

