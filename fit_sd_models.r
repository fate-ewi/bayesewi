
library(dplyr)
library(bayesewi)

# This block loads in the data, and fits EWI models to each
da_data = read.csv("/users/eric.ward/downloads/da.examples.csv")

da_data = select(da_data, example, dec.date, DA.ppm) %>%
  rename(y = DA.ppm, x = dec.date) %>%
  filter(y != 0)
da_data$y = log(da_data$y)

models = list()
for(i in unique(da_data$example)) {
  df = filter(da_data, example==i) %>%
    select(-example)
  df$x = round(365*(df$x - min(df$x)))
  models[[i]] = fit_ewi(df, ewi_model="sv", iter = 2000, chains=3)
}
saveRDS(models, file="fitted_sv_models.rds")

models = readRDS(file="fitted_sv_models.rds")
m = models[[1]]

pdf("DA_ewi.pdf")
for(i in 1:8) {
  plot_estimates(models[[i]], alpha = 0.05, link_space=TRUE)
}
dev.off()


