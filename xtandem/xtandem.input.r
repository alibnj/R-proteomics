# This code fixes the paths in the .xml file using regex. The .xml files are from X! TANDEM

#------------- FOR HUMAN -------------------------------------------

for (i in 1:129)
{
  mzxml <- paste("/gss_gpfs_scratch/ali.b/xtandem/single.cells.v2/human/h", i, ".mzXML", sep = "")
  result <- paste("/gss_gpfs_scratch/ali.b/xtandem/single.cells.v2/human/", i, ".tandem.xml", sep = "")
  logg <- paste("/gss_gpfs_scratch/ali.b/xtandem/single.cells.v2/human/log", i, ".txt", sep = "")
  input <- paste("C:/Users/Ali/Desktop/xtandem/v2/human/input", i, ".xml", sep = "")
  
  tx  <- readLines("C:/Users/Ali/Desktop/xtandem/v2/human/input.human.xml")
  
  tx2 <- gsub(pattern = "/path_to_mzxml/", replace = mzxml, x = tx)
  tx3  <- gsub(pattern = "/path_to_output/", replace = result, x = tx2, fixed=TRUE)
  tx4  <- gsub(pattern = "/path_to_log/", replace = logg, x = tx3, fixed=TRUE)
  
  writeLines(tx4, con=input)
}

#------------- FOR MOUSE -------------------------------------------

for (i in 1:15)
{
  mzxml <- paste("/gss_gpfs_scratch/ali.b/xtandem/single.cells.v2/mouse/m", i, ".mzXML", sep = "")
  result <- paste("/gss_gpfs_scratch/ali.b/xtandem/single.cells.v2/mouse/", i, ".tandem.xml", sep = "")
  logg <- paste("/gss_gpfs_scratch/ali.b/xtandem/single.cells.v2/mouse/log", i, ".txt", sep = "")
  input <- paste("C:/Users/Ali/Desktop/xtandem/v2/mouse/input", i, ".xml", sep = "")
  
  tx  <- readLines("C:/Users/Ali/Desktop/xtandem/v2/mouse/input.mouse.xml")
  
  tx2 <- gsub(pattern = "/path_to_mzxml/", replace = mzxml, x = tx)
  tx3  <- gsub(pattern = "/path_to_output/", replace = result, x = tx2, fixed=TRUE)
  tx4  <- gsub(pattern = "/path_to_log/", replace = logg, x = tx3, fixed=TRUE)
  
  writeLines(tx4, con=input)
}
