#Creadora: Michelle Clempner
#Gráfica de Volcano



library(pacman)

p_load("readr",     #Llamar las bases de datos
       "ggplot2",   #Graficar 
       "ggrepel",   #Etiquetar datos en una gráfica
       "matrixTests", #Realizar la prueba estadística
       "dplyr")      # Facilita el manejo de datos

datos <- read_csv(file= "https://raw.githubusercontent.com/ManuelLaraMVZ/Transcript-mica/main/datos_miRNAs.csv")
head(datos)

# Extracción de genes controles (referencia)

controles <- datos %>% 
  filter(Condicion=="Control")

head(controles)

# Sacar los promedios 

promedios_controles <- controles %>% 
  summarise(Mean_C1=mean(Cx1),
            Mean_C2=mean(Cx2),
            Mean_C3=mean(Cx3),
            Mean_T1=mean(T1),
            Mean_T2=mean(T2),
            Mean_T3=mean(T3)) %>% 
  mutate(Gen="promedio_controles") %>% 
  select(7,1,2,3,4,5,6)
promedios_controles            

#Extraer los genes de la tabla "datos"

genes <- datos %>% 
  filter(Condicion=="Target") %>% 
  select(-2)
head(genes)


#Sacar el 2^-DCT

DCT <- genes %>% 
  mutate (DCT_C1=2^-(Cx1-promedios_controles$Mean_C1),
          DCT_C2=2^-(Cx2-promedios_controles$Mean_C2),
          DCT_C3=2^-(Cx3-promedios_controles$Mean_C3), 
          DCT_T1=2^-(T1-promedios_controles$Mean_T1),
          DCT_T2=2^-(T2-promedios_controles$Mean_T2),
          DCT_T3=2^-(T3-promedios_controles$Mean_T3)) %>% 
  select(-2,-3,-4,-5,-6,-7)


DCT

promedio_genes <- DCT %>% 
  mutate(Mean_DCT_Cx=(DCT_C1+DCT_C2+DCT_C3)/3,
         Mean_DCT_Tx=(DCT_T1+DCT_T2+DCT_T3)/3)  


promedio_genes


#Definir prueba estadística: t-Student

tvalue_gen <- row_t_welch(promedio_genes[,c("DCT_C1",
                                            "DCT_C2",
                                            "DCT_C3")],
                          promedio_genes[,c("DCT_T1",
                                            "DCT_T2",
                                            "DCT_T3")])
View(tvalue_gen)

FCyPV <- promedio_genes %>% 
  select(1,8,9) %>% 
  mutate(p_value = tvalue_gen$pvalue,
         Fold_change =Mean_DCT_Tx/Mean_DCT_Cx) %>% 
  select(1,4,5)

FCyPV

Logs <- FCyPV %>% 
  mutate(LPV = -log10(p_value),
         LFC = log2(Fold_change)) %>% 
  select (1,4,5)

Logs


vulcano <- ggplot(Logs, mapping = aes(x=LFC,
                                      y=LPV))+
  geom_point(size=2,
             color="gray")+
  theme_classic()+
  labs(title = "Análisis comparativo de miRNAs",
       caption = "Creadora: Michelle Clempner",
       x= "Log2 (2-DDCT)",
       y="log10(valor de p)")

vulcano



#Límites

limite_p <- 0.05

limite_FC <- 1.5

down_regulated <- Logs %>% 
  filter(LFC < -log2(limite_FC),
         LPV > -log10(limite_p))

down_regulated

up_regulated <- Logs %>% 
  filter(LFC > log2(limite_FC),
         LPV > -log10(limite_p))
up_regulated

top_down_regulated <- down_regulated %>% 
  arrange(desc(LPV)) %>% 
  head(5)

top_down_regulated

top_up_regulated <- up_regulated %>% 
  arrange(desc(LPV)) %>% 
  head(5)

top_up_regulated


#Mejorar la Gráfica

vulcano2 <- vulcano+
  geom_hline(yintercept = -log10(limite_p),
             linetype = "dashed") +
  geom_vline(xintercept = c(-log2(limite_FC), log2 (limite_FC)),
             linetype = "dashed")

vulcano2

vulcano3 <- vulcano2+
  geom_point(data= up_regulated,
             x= up_regulated$LFC,
             y= up_regulated$LPV,
             alpha=1,
             size=3,
             color= "#E64B35B2" )+
  geom_point(data= down_regulated,
             x= down_regulated$LFC,
             y= down_regulated$LPV,
             alpha=1,
             size=3,
             color= "#3C5488B2" )
vulcano3

vulcano4 <- vulcano3 +
  geom_label_repel(data= top_up_regulated,
                   mapping = aes(x= top_up_regulated$LFC,
                                 y= top_up_regulated$LPV),
                   label= top_up_regulated$Gen,
                   max.overlaps = 50,
                   box.padding = 0.5, # Ajusta estos valores según tu preferencia
                   point.padding = 1, # Ajusta estos valores según tu preferencia
                   segment.curvature = 0.2) + # Ajusta estos valores según tu preferencia
  geom_label_repel(data= top_down_regulated,
                   mapping = aes(x= top_down_regulated$LFC,
                                 y= top_down_regulated$LPV),
                   label= top_down_regulated$Gen,
                   max.overlaps = 50,
                   box.padding = 0.5, # Ajusta estos valores según tu preferencia
                   point.padding = 1, # Ajusta estos valores según tu preferencia
                   segment.curvature = 0.2) # Ajusta estos valores según tu preferencia

vulcano4  

ggsave("Vulcano4.jpeg",
       plot = vulcano4,
       height = 5, width =6,
       dpi = 300)

