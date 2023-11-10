# M <- file.path(Sys.getenv("HOME"),".R", ifelse(.Platform$OS.type == "windows", "Makevars.win", "Makevars"))
# file.edit(M)
# system('wmic cpu get caption, deviceid, name, numberofcores')

dotR <- file.path(Sys.getenv("HOME"), ".R")
if (!file.exists(dotR)) dir.create(dotR)
M <- file.path(dotR, "Makevars.win")
if (!file.exists(M)) file.create(M)
cat("\nCXX14FLAGS=-O3 -march=corei7 -mtune=corei7",
    "CXX14 = $(BINPREF)g++ -m$(WIN) -std=c++1y",
    "CXX11FLAGS=-O3 -march=corei7 -mtune=corei7",
    file = M, sep = "\n", append = TRUE)

file.edit(M)
