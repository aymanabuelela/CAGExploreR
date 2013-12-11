html.report <-
function(x,k,fig.dir = "Figures",report.name="Switch Report")
{

dir.create(fig.dir)
oldwd = getwd()
#setwd(fig.dir)
  .HTML.file <<- file.path(getwd(),paste(report.name,".html",sep=""))
  HTML(as.title("Promoter Switching Report"),append=FALSE)
setwd(fig.dir)
  HTMLhr()
	counter = 1
for(GENE in k)
{
	png(paste(as.character(GENE),".png",sep=""),width=1500,height=750,units="px",pointsize=20)
	PROMS = plotcomp(x,GENE)
	dev.off()
	HTML(paste(counter,"HGNC:",GENE,", promoters:",sep=" "))
	HTML(PROMS)
	HTMLInsertGraph(paste(getwd(),"/",GENE,".png",sep=""),WidthHTML=1500,HeightHTML=750)
	HTMLhr()
	counter = counter + 1
}
 setwd(oldwd)

}
