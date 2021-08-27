# cargar libs
library("dplyr")
library("ggplot2")
library("ggsci")

## Visualizar Rango Dinamico ====
# creamos un datframe para el rango dinamico 1
rango1 <- final_data_b %>% 
  select( id.cond1, log10.cond1 ) %>% 
  rename( IDs = id.cond1,
          log10.average = log10.cond1 ) %>% 
  mutate( condition = condicion1 )

# creamos un datframe para el rango dinamico 2
rango2 <- final_data_b %>% 
  select( id.cond2, log10.cond2 ) %>% 
  rename( IDs = id.cond2,
          log10.average = log10.cond2 ) %>% 
  mutate( condition = condicion2 )

# Juntamos ambos rangos
rango_final <- bind_rows( rango1, rango2 )

# graficamos con puntos
dinamico1 <- ggplot( data = rango_final,
                     mapping = aes(x = IDs,
                                   y = log10.average,
                                   color = condition ) ) +
  geom_point( alpha = 0.5 )

# Visualizamos
dinamico1

# cambiamos colores y titulos
dinamico2 <- dinamico1 +
  scale_color_d3() +
  labs( title = "Dynamic Range",
        x = "IDs",
        y = "log10(averages)")

# Visualizamos
dinamico2

# mejoramos el tema
dinamico3 <- dinamico2 +
  theme_classic( base_size = 15 ) +
  theme( plot.background = element_rect( fill = "white" ))

# Visualizamos
dinamico3

# Guardamos el rango
# guardamos plot
guardar_imagen( nombre_imagen = "rango_dinamico.png",
                el_plot = dinamico3, ancho = 10, alto = 5 )

# FIN DEL SCRIPT