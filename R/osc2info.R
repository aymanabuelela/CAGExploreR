osc2info <-
function(x) #convert OSC prom names to component information
{
  split1 = strsplit(x,",",fixed=T)
  strand = as.character(unlist(lapply(split1,function(x) x[2])))
  temp.anot = unlist(lapply(split1,function(x) x[1]))
  split1 = strsplit(temp.anot,":",fixed=T)
  chr = as.character(unlist(lapply(split1,function(x) x[1])))
  temp.anot = unlist(lapply(split1,function(x) x[2]))
  split1 = strsplit(temp.anot,"..",fixed=T)
  start = as.numeric(unlist(lapply(split1,function(x) x[1])))
  end = as.numeric(unlist(lapply(split1,function(x) x[2])))
  return(list(chr=chr,strand=strand,start=start,end=end))
}
