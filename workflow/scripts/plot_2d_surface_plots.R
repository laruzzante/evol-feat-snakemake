df <- merged_orthogroup_features

library(ggplot2)

ggplot(df, aes(x=UNI, y=RUN) ) +
  geom_bin2d(bins=70) +
  scale_fill_continuous(type = "viridis") +
  theme_bw()


library(plot3D)

x <- df$AGE
y <- df$RUN

##  Create cuts:
n_cuts <- 30
x_c <- cut(x, n_cuts)
y_c <- cut(y, n_cuts)

##  Calculate joint counts at cut levels:
z <- table(x_c, y_c)

##  Plot as a 3D histogram:
hist3D(z=z, border="black", xlab='AGE', ylab='RUN', zlab='counts')

##  Plot as a 2D heatmap:
image2D(z=z, border="black")
