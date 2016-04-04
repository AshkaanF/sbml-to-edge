library(SBMLR)

sbml.to.edge.list <- function(sbml.object){
  x <- sbml.object
  
  reactions.list <- list()
  
  ## loop through each reaction in the metabolic model
  for(i in 1:length(x$reactions)){
    temp <- x$reactions[i]
    reactants <- unlist(temp[[1]][3])
    products <- unlist(temp[[1]][4])
    
    ## store results as list element
    temp.edge.list <- as.data.frame(matrix(ncol = 4, nrow = length(reactants) * length(products)))
    names(temp.edge.list) <- c('reactant', 'product', 're.comp', 'prod.comp')
    
    ## add edges
    temp.edge.list$reactant <- rep(reactants, each = length(products))
    temp.edge.list$product <- rep(products)
    
    reactions.list[[i]] <- temp.edge.list
    
    ## drop us a line
    if(i == 1){bar.vec <- c(na.omit(seq(1:length(x$reactions))[1:length(x$reactions) * round(length(x$reactions) / 10)]))
               cat('|')}
    if(i %in% bar.vec == TRUE){cat('=====|')}
  }
  
  ## print edge list after binding reactions.list
  final <- do.call('rbind', reactions.list)
  final <- final[-which(final$product == '0'), ]
  final <- final[-which(final$product == '1000'), ]
  
  cat('\n...wrapping things up\n')
  
  ## reformat node labels
  for(q in 1:length(final$reactant)){
    final$re.comp[q] <- unlist(strsplit(final$reactant[q], '_'))[2]
    final$prod.comp[q] <- unlist(strsplit(final$product[q], '_'))[2]
  }
  
  ## subset by external compartment and remove biomass component
  final <- subset(final, re.comp == 'c0')
  final <- subset(final, prod.comp != 'b')
  unique(final)
}