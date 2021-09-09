# cargar librerias
library( "factoextra" )
library( "dplyr" )
library( "stringr" )
library( "scales" )
library( "ggsci" )
library( "tidyr" )
library( "pheatmap" )

# desactivamos toda la notacion cientifica
options( scipen = 666 )

# juntamos las protes up y down
base_pca <- bind_rows( proteinas_up, proteinas_down )

# ponemos los Acc como nombres de filas
rownames( base_pca ) <- base_pca$Accession

# volvemos a sacar lo siguiente: 
# Vamos a buscar el numero de columnas que corresponden a las muestras de la condicion 1
re_posicion_muestras_1 <- str_detect( string = colnames( base_pca ),
                                      pattern = expresion_1 ) %>%
  which()

# vamos a extraer un dataframe que contenga solo las columnas de las inyecciones
re_data_condicion_1 <- base_pca %>% 
  select( re_posicion_muestras_1 )

# volvemos a sacar lo siguiente: 
# Vamos a buscar el numero de columnas que corresponden a las muestras de la condicion 1
re_posicion_muestras_2 <- str_detect( string = colnames( base_pca ),
                                      pattern = expresion_2 ) %>%
  which()

# vamos a extraer un dataframe que contenga solo las columnas de las inyecciones
re_data_condicion_2 <- base_pca %>% 
  select( re_posicion_muestras_2 )

# Juntamos los datas que contienen SOLO las columnas de medicion de inyecciones
base2_pca <- bind_cols( re_data_condicion_1, re_data_condicion_2 )

# Hacemos un transpose de columnas a filas
trans_pca <- t( base2_pca ) %>% 
  as.data.frame( )

# agregamos columna de nombre de muestra limpio, y condicion a la que pertenece
trans2_pca <- trans_pca %>% 
  mutate( muestra = row.names(.) %>% 
            str_sub( start = 2 ) ) %>% 
  left_join( x = .,
             y = muestras,
             by = "muestra" )

# el join pierde los nombres, asi que los recuperamos
rownames( trans2_pca ) <- trans2_pca$muestra

# comenzamos con el pca
# cual es el maximo de columnas en trans 2? - Esto camibara dependiendo del numero de proteinas que entran al analisis
ultima_col <- ncol(trans2_pca)
penultima_col <- ultima_col - 1

# haremos el pca con trans2 pero sin sus ultimas dos columnas
pca_resultados <- trans2_pca %>%
  select( -ultima_col, -penultima_col ) %>% 
  prcomp( scale = TRUE )

# Visualizamos el screeplot
maxy = 100
marcas_y = seq( from = 0, to = maxy, by = 10 )
etiquetas_y = (marcas_y / 100) %>% percent( )

screeplot <- fviz_eig( pca_resultados ) +
  scale_y_continuous( breaks = marcas_y,
                      labels = etiquetas_y ) +
  theme_classic( base_size = 15 )

# Vis
screeplot

# guardar el scree
guardar_imagen( nombre_imagen = "screeplot.png",
                el_plot = screeplot, ancho = 7, alto = 5 )

# Scamos el primer PCA
pca_ind <- fviz_pca_ind( pca_resultados,
                         axes = c(pc_enx, pc_eny),
                         geom = "point", pointsize = 3,
                         col.ind = as.factor( trans2_pca[, ultima_col] ) ,
                         invisible = "quali" ) +
  scale_color_d3() +
  theme_bw( base_size = 15 ) +
  theme( legend.title = element_blank() )

# vis
pca_ind

# guardar el pcaind
guardar_imagen( nombre_imagen = "pca.png",
                el_plot = pca_ind, ancho = 7, alto = 5 )

# graficamos el biplot
bip <- fviz_pca_biplot( pca_resultados,
                        axes = c(pc_enx, pc_eny),
                        geom.ind = "point", pointsize = 3,
                        geom.var = "arrow", col.var = "black", alpha.var = 0.3,
                        col.ind = as.factor( trans2_pca[, ultima_col] ) ,
                        invisible = "quali" ) +
  scale_color_d3() +
  theme_bw( base_size = 15 ) +
  theme( legend.title = element_blank() )

# vis
bip

# guardar el biplot
guardar_imagen( nombre_imagen = "biplot.png",
                el_plot = bip, ancho = 7, alto = 5 )

# Vamos a generar un PCP (parallel coordinate plot)
# Permite ver varios ejes del PCA al mismo tiempo
# Sacamos los resultados completos
# Results for individuals
resultados_individuales <- get_pca_ind( pca_resultados )

todos_pc <- resultados_individuales$coord %>% 
  as.data.frame( ) %>% 
  mutate( muestra = row.names(.) ) %>% 
  left_join( x = .,
             y = muestras,
             by = "muestra" )

# volvemos a calcular la ultima columna
re_ultima_col <- ncol( todos_pc )
re_penultima <- re_ultima_col - 1

# pasamos los PC a formato largo
pc_largo <- todos_pc %>% 
  pivot_longer( cols = -c(re_penultima:re_ultima_col),
                names_to = "PC",
                values_to = "coordinate" )

# Hagamos el PCP con lineas y puntos
pcp <- ggplot( data = pc_largo,
               mapping = aes( x = PC,
                              y = coordinate,
                              color = condition,
                              fill = condition) ) +
  geom_line( mapping = aes( group = muestra ), size = 1 ) +
  geom_point( shape = 21, size = 3, color = "black" ) +
  labs( title = "Parallel Coordinate Plot" ) +
  theme_light( base_size = 15 ) +
  theme( panel.grid.minor.y = element_blank(),
         legend.position = "top",
         legend.title = element_blank() )

# vis
pcp

guardar_imagen( nombre_imagen = "PCP.png",
                el_plot = pcp, ancho = 7, alto = 5 )

## HEATMAP ====
# convertimos a matriz
heat_data <- as.matrix( base2_pca  )

# escalamos los datos
heat_escalado <- scale( heat_data )

# hacemos el pheat
# vamos a etiqutear columnas y filas...
# Cargar el etiquetado de columnas
col_ids <- muestras               # Solo dice el nombre de la columna, a que estado celular pertenece
rownames( col_ids ) <- paste0("X", col_ids$muestra)
# dejar solo los nombres de etiquetas
# Lo requiere pheatmap
etiquetas_col <- col_ids %>%
  select( condition )

# etiquetamos protes
# Lo utilizara para ponerle color al gen
etiquetas_prote <- bind_rows( proteinas_up %>% mutate( state = "UP"),
                              proteinas_down %>% mutate( state = "DOWN") )

# ponemos los Acc como nombres de filas
rownames( etiquetas_prote ) <- etiquetas_prote$Accession

# solo nos quedamos con la etiqueta
etiquetas_prote2 <- etiquetas_prote %>% 
  select( state )

# Se deben definir varias escalas de colores al mismo tiempo
# lo guarda en una lista
# que a su vez contiene dos vectores
mis_colores = list( condition = c( "#1f77b4", "#ff7f0e"),
                    state = c("UP" = "#d9352a", "DOWN" = "#4979b6") )  

# aqui corregimos un bug donde no se puede pasar un objeto como nombre del elemento en vector condition
names( mis_colores[[1]] ) <- c(condicion1, condicion2)

heat1 <- pheatmap( mat = heat_escalado,
                   main = subtitulo,
                   scale = "row",
                   fontsize_row = 5,
                   fontsize_col = 5,
                   cutree_rows = 2,
                   cutree_cols = 2,
                   show_rownames = TRUE,
                   annotation_colors = mis_colores,
                   annotation_row = etiquetas_prote2,
                   annotation_col = etiquetas_col )

# guardar el heatmap
guardar_imagen( nombre_imagen = "heatmap.png",
                el_plot = heat1, ancho = 10, alto = 10 )

# FIN DEL SCRIPT
