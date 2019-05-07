@echo off

Rscript R\beta_bias.R 3 0.5
Rscript R\Catch_bias.R 3 0.5
Rscript R\H_bias.R 3 0.5

Rscript R\beta_bias.R 4 0.5
Rscript R\Catch_bias.R 4 0.5
Rscript R\H_bias.R 4 0.5