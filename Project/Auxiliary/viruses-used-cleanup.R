###
# Cleaning up "viruses_used.Rds"
# Zane Billings
# 2024-04-03
# Literally just writes the viruses_used.Rds file from Amanda's project to
# a CSV so I can open it in excel.
###

box::use(
	readr,
	here
)

dat <- readr::read_rds(here::here("data/viruses_used.rds"))

readr::write_csv(dat, here::here("Auxiliary/viruses-used.csv"))

# END OF FILE ####
