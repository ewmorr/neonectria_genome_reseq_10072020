require(prism)

prism_set_dl_dir("/mnt/home/garnas/ericm/PRISM_dailies_ppt")
get_prism_dailys(type = "ppt", minDate="2007-01-01", maxDate="2018-12-31", keepZip = F)

