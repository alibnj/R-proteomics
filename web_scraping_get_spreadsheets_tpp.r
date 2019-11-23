install.packages('RSelenium')
library('RSelenium')


for (i in 101:129) {
  link <- paste("http://localhost:10401/tpp/cgi-bin/PepXMLViewer.cgi?page=1&columns=Pprobability%2CGspectrum%2CSexpect%2CGions%2CGpeptide%2CGprotein%2CGcalc_neutral_pep_mass%2CGretention_time_sec%2CGindex%2CGstart_scan%2CGassumed_charge%2CGprecursor_neutral_mass%2CGMZratio%2CGpI%2CGcompensation_voltage%2CGprecursor_intensity%2CSxcorr%2CSdeltacn%2CSdeltacnstar%2CSspscore%2CSsprank%2CGions_old%2CGnum_tol_term%2CGnum_missed_cleavages%2CGmassdiff%2CPfval%2CQlibra1%2CQlibra2%2CQlibra3%2CQlibra4%2CQlibra5%2CQlibra6%2CQlibra7%2CQlibra8%2CQlibra9%2CQlibra10&displayState=otherActionsDiv&exportSpreadsheet=1&sortField2=Gspectrum&sortDir2=0&FmPprobability2=&FMPprobability2=&xmlFileName=c%3A%2FInetpub%2Fwwwroot%2FISB%2Fdata%2Fcomet-p07-libra%2Fh", i, ".P0.7.xml&perPage=50&sortField=Gspectrum&sortDir=0&highlightedPeptideText=&highlightedProteinText=&highlightedSpectrumText=&libraResultType=absolute&expandProteinList=condense&minimizeTableHeaders=yes&requiredAA=&requiredPeptideText=&requiredProteinText=&requiredSpectrumText=&FmGnum_tol_term=&FMGnum_tol_term=&FmGnum_missed_cleavages=&FMGnum_missed_cleavages=&fm1Sxcorr=&fm1Sdeltacn=&fM1Ssprank=&fm2Sxcorr=&fm2Sdeltacn=&fM2Ssprank=&fm3Sxcorr=&fm3Sdeltacn=&fM3Ssprank=&FmSexpect=&FMSexpect=&FmPprobability=&FMPprobability=&jumpPage=1", sep = "")
  shell.exec(link) #Opens the url using the default browser 
  }



for (j in 20:32)
{
  #checkForServer()
  #startServer()
  remDrv <- remoteDriver()
  remDrv$open()
  C <- data.frame(Gene = character(), Organism = character(), Locus = character(), Region = character(), Length = numeric())
  filename <- paste("C:/Users/banijamali.s/Dropbox/Genetics/R/New Results/UTRdb Data Extraction/185.8/UTRdb.185.8.rem", j, ".csv", sep = " ")
  for (i in ((j*1500)+1):((j+1)*1500))
  {
    remDrv$navigate('http://utrdb.ba.itb.cnr.it/search')
    remDrv$findElement(using = "xpath", "//select[@id = 'utr_db']/option[@value = '1']")$clickElement()
    remDrv$findElement(using = "xpath", "//input[@id ='ids']")$sendKeysToElement(list(m.185.8.rem[i,1], "\uE007"))
    doc <- htmlParse(remDrv$getPageSource()[[1]])
    tab <- readHTMLTable(doc)
    tab <- as.data.frame(tab)
    C <- rbind(C,tab)
    Sys.sleep(1)
  }
  write.csv(C, filename)
  remDrv$closeWindow()
  remDrv$quit()
  remDrv$closeServer()
  Sys.sleep(5)
}
