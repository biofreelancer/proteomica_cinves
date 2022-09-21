# cargar librerias
library("dplyr")
library("ggplot2")
library("cowplot")
library("ggsci")
library("ggrepel")

## Visualizar Volcano ====
# filtramos la data base del volcano
base_volcano <- final_data_b %>% 
  filter( Anova..p. < pvalue_limite,
          Peptide.count >= peptide_count_limite,
          Unique.peptides >= unique_peptide_limite,
          Count.cond1 >= count_limite,
          Count.cond2 >= count_limite,
          CV.cond1 > cv_limite,
          CV.cond2 > cv_limite )

# Dibujamos el volcano base con todos los puntos
volcano1 <- ggplot( data = final_data_b,
                    mapping = aes( x = log2.Ratio,
                                   y = nlog10.Anova ) ) +
  geom_point( alpha = 0.3,                                      # Usamos una transparencia del 10%
              shape = 1,                                        # la forma es de circulo hueco
              color = "#2E294E" )

# Visualizamos
volcano1

# Dibujamos las lineas en los limites de pvalue y de log2
vertical_lines.v <- c(-ratio_limite, ratio_limite)
  
volcano2 <- volcano1 +
  geom_hline( yintercept = -log10( pvalue_limite ),
              linetype="dashed" ) +
  geom_vline( xintercept = vertical_lines.v,
              linetype="dashed")

# Visualizamos
volcano2

# Pintamos la capa para las proteinas UP
# Sacamos las proteinas up 
proteinas_up <- base_volcano %>% 
  filter( log2.Ratio > ratio_limite )

# graficamos las up
volcano3 <- volcano2 +
  geom_point( data = proteinas_up,
              color = "#2D728F",
              alpha = 0.7,
              size = 3 )

# Visualizamos
volcano3

# Pintamos la capa para las proteinas DOWN
# Sacamos las proteinas down
proteinas_down <- base_volcano %>% 
  filter( log2.Ratio < -ratio_limite )

# graficamos las down
volcano4 <- volcano3 +
  geom_point( data = proteinas_down,
              color = "#AB3428",
              alpha = 0.7,
              size = 3 )

# Visualizamos
volcano4

# Agregamos mas dientes (breaks) al eje X
# definimos el valor maximo absoluto que encontramos log2 ratio
absolute_max_x <- final_data_b$log2.Ratio %>%
  abs( ) %>%
  max() %>%
  ceiling()

# creamos valores para el eje x
limites <- c( -absolute_max_x, absolute_max_x)

# marcas eje X
ejex <- data.frame( marcas = -absolute_max_x:absolute_max_x) %>% 
  mutate( fila = 1:nrow( . )  ) %>% 
  mutate( espar = ifelse( test = fila %% 2 == 0,
                          yes = "par",
                          no = "non" ) ) %>% 
  mutate( etiquetas = ifelse( test = espar == "par",
                              yes = marcas,
                              no = "" ) )

# graficamos
volcano5 <- volcano4 +
  scale_x_continuous( limits = limites,
                      breaks = ejex$marcas,
                      labels = ejex$etiquetas )

# Visualizamos
volcano5

# creamos titulos, subs, etc
subtitulo <- paste( condicion2, "vs", condicion1)

titulo_eje_x <- paste( "log2(",
                       condicion2,
                       "/",
                       condicion1,
                       ")" )

# agregamos titulos, subs, etc
volcano6 <- volcano5 +
  labs(title = "Peptidoma diferencial",
       subtitle = subtitulo,
       x = titulo_eje_x ,
       y = "-log10( ANOVA p-value )" )

# Visualizamos
volcano6

# giramos coordenadas, limpiamos tema
volcano7 <- volcano6  +
  coord_flip( ) +
  theme_half_open( font_size = 15 ) +
  theme( plot.background = element_rect( fill = "white" ) )

# Visualizamos
volcano7

# guardamos plot
guardar_imagen( nombre_imagen = "volcano.png",
                el_plot = volcano7, ancho = 10, alto = 5 )

# # sacar los N mas up
top_up <- proteinas_up %>%
  arrange( log2.Ratio ) %>%
  tail( n = n_best )

# sacar los N mas down
top_down <- proteinas_down %>%
  arrange( log2.Ratio ) %>%
  head( n = n_best )

# agregamos las etiquetas al plot
volcano8 <- volcano7 +
geom_label_repel( data = top_up,
                  mapping = aes( label = Accession ),
                  max.overlaps = 50 ) +
geom_label_repel( data = top_down,
                  mapping = aes( label = Accession ),
                  max.overlaps = 50 )

# Visualizamos
volcano8

# Guardamos el plot
# guardamos plot
guardar_imagen( nombre_imagen = "volcano_etiquetado.png",
                el_plot = volcano8, ancho = 10, alto = 5 )

# FIN DE SCRIPT
