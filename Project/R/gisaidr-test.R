###
#  Testing GISAIDR since GISAID won't let me look at the flu sequences online
# Zane
# 2024-03-26
###

library(GISAIDR)

username <- Sys.getenv("GISAIDR_USERNAME")
password <- Sys.getenv("GISAIDR_PASSWORD")

credentials <- login(username = username, password = password, database = "EpiFlu")
