
###------------------------------------------------------------------##
###       Codes for installing cmdstanR
###       If you're using R 4.x.x please see the discussion at 
###       [https://discourse.mc-stan.org/t/cmdstanr-installation-error-command-not-found-win-processx-c-977/11970/4]
###       Also, my question and Rok's answer at:
###       [https://discourse.mc-stan.org/t/cmdstanr-error-cc1plus-exe-sorry-unimplemented-64-bit-mode-not-compiled-in/16720/6]
###------------------------------------------------------------------##


# install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
library(cmdstanr)
library(posterior)
library(bayesplot)
color_scheme_set("brightblue")

## Install Rtools, set path correctly, then run this: 

check_cmdstan_toolchain(fix = TRUE)

install_cmdstan(overwrite = TRUE)
# cmdstan_make_local(cpp_options = "CXXFLAGS += --Wno-int-in-bool-context", append = TRUE)

file <- file.path(cmdstan_path(), "examples", "bernoulli", "bernoulli.stan")
mod <- cmdstan_model(file)

# install.packages(
#   "https://win-builder.r-project.org/gdaEAY9p8V7I/StanHeaders_2.21.0-6.zip",
#   repos = NULL, type = "win.binary")
# install.packages(
#   "https://win-builder.r-project.org/iT79L894j228/rstan_2.21.2.zip",
#   repos = NULL, type = "win.binary")

data_list <- list(N = 10, y = c(0,1,0,0,0,0,0,0,0,1))

fit <- mod$sample(
  data = data_list,
  seed = 123,
  chains = 4,
  parallel_chains = 2,
  refresh = 500
)

fit$summary()
# pkgbuild::has_build_tools(debug = TRUE)

