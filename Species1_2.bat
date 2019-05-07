@echo off

Rscript R\beta_bias.R 1 0.5
Rscript R\Catch_bias.R 1 0.5
Rscript R\H_bias.R 1 0.5

Rscript R\beta_bias.R 2 0.5
Rscript R\Catch_bias.R 2 0.5
Rscript R\H_bias.R 2 0.5


