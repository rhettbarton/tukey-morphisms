#install.packages(c("igraph","R.utils","gtools","data.table","RcppGreedySetCover"))
#install.packages("BiocManager")
#BiocManager::install("Rgraphviz")

library(igraph)
library(R.utils)
library(gtools)
library(data.table)
library(RcppGreedySetCover)
library(sqldf)
#library(hasseDiagram)

#Set Working Directory to save the data and output
setwd("/bsuhome/rbarton/")

#Save Workspace
save.image(file = "TukeyMorphism_start.RData")

#####################################################################
################# Begin Function Definition #########################
#####################################################################

graph_from_matrix <- function(M,plotGraph = T,graphName = "",saveGraph = F) {
  #Convert incidence matrices to graphs
  library(igraph)
  graphs_temp <- as.directed(graph.incidence(M), mode = "arbitrary")
  if(plotGraph == T) {
    if(length(M) == 1) {
      i <- 1
      j <- 1
    } else {
      i <- nrow(M)
      j <- ncol(M)
    }
    layout <- matrix(c(seq(0,i-1),seq(0,j-1),rep(0,i),rep(1,j)),byrow = F,nrow = i+j, ncol = 2)  
    plot(graphs_temp,layout = layout, vertex.color = 'black', vertex.label = NA, vertex.size = 7, 
         edge.arrow.size=1,main = paste(graphName,sep=''))
  }
  if(saveGraph == T) {
    assign(graphName,graphs_temp,envir = .GlobalEnv)
  }
  return(graphs_temp)
}

generate_finite_relations <- function(Usizes,Vsizes,
                                      saveMatrices = T,saveGraphs = F,
                                      matricesListName = 'matrices',graphsListName = 'graphs'){
  library(R.utils)
  #Generate finite relations of a given size
  library(igraph)
  library(gtools)
  grid <- expand.grid(Usizes,Vsizes)
  #Create all possible binary vectors of length j
  #2^j possibilities at each step
  perm <- lapply(X = Vsizes,function(k) permutations(n= 2, r = k, v = 0:1, repeats.allowed = TRUE)) 
  #Get "all" ixj matrices by combining (with replacement) the vectors from the previous step i times
  #The zero vectors can be left out, since we are only interested in relations that have dominating families 
  #((2^j-1)+i-1)!/(i!((2^j-1)-1)!) possibilities at each step    
  combi <- lapply(X = 1:nrow(grid),function(k) combinations(n = nrow(perm[[match(grid[k,2],Vsizes)]])-1, r = grid[k,1], 
                                                            v = 2:nrow(perm[[match(grid[k,2],Vsizes)]]), repeats.allowed = TRUE))
  matrix_temp <- lapply(X = 1:length(combi),
                        function(k) lapply(X = 1:nrow(combi[[k]]), 
                                           function(n) matrix(perm[[match(grid[k,2],Vsizes)]][combi[[k]][n,],],
                                                              byrow = F, nrow = grid[k,1], ncol = grid[k,2])))
  if(saveMatrices == T) {
    names(matrix_temp) <- lapply(X = 1:nrow(grid), function(k) paste0("Size_",grid[k,1],"x",grid[k,2]))
    matrix_name <- paste(matricesListName)
    assign(matrix_name,matrix_temp, envir = .GlobalEnv)
  }
  if(saveGraphs == T) {
    graphs_temp <- lapply(X = matrix_temp,function(M) as.directed(graph.incidence(M), mode = "arbitrary"))        
    graphs_name <- paste(graphsListName)
    assign(graphs_name,graphs_temp, envir = .GlobalEnv)
  }      
}  

count_of_graph_generation <- function(Usizes,Vsizes) {
  n <- 0
  for(j in Vsizes){
    #2^j possibilities at each step
    for(i in Usizes) {
      #Get "all" ixj matrices by combining (with replacement) the vectors from the previous step i times
      #The zero vectors can be left out, since we are only interested in relations that have dominating families 
      #((2^j-1)+i-1)!/(i!((2^j-1)-1)!) possibilities at each step 
      n <- n + factorial((2^j-1)+i-1)/(factorial(i)*factorial((2^j-1)-1))
    }
  }
  return(n)  
}

dual_relation <- function(A) {
  #Find the dual of a finite relation
  Adual <- t(1- A)
  return(Adual)
}

setsystem_from_matrix <- function(M){
  library(data.table)
  #Convert incidence matrix to two column data frame (set system): 
  j = ncol(M)
  S1 <- c()
  S2 <- c()
  for (k in 1:j) {
    S1 <- c(S1,rep(k,length(which(M[,k] == 1)))) # - 1st column is points in A+ (columns in the incidence matrix) 
    S2 <- c(S2,which(M[,k] == 1)) # - 2nd column contain the points in A- that they relate to     
  }
  S <- data.table("set" = S1, "element" = S2)
  return(S)
}

dominating_number <- function(M) {
  #Find the (approximate) dominating number of a finite relation 
  #Test for if a dominating number exists
  if(0 %in% rowSums(M)) {return(NA)}
  #Test for dominating number = 1
  if(nrow(M) %in% colSums(M)) {return(1)}
  #Test for dominating number = 2
  combs <- t(combn(x = 1:ncol(M),m = 2))
  comp_value <- apply(X = combs, MARGIN = 1, function(ind) pmax(M[,ind[1]],M[,ind[2]]))
  nrow(M) %in% colSums(comp_value)
  if(nrow(M) %in% colSums(comp_value)) {return(2)}
  #Test for dominating number = 3
  combs <- t(combn(x = 1:ncol(M),m = 3))
  comp_value <- apply(X = combs, MARGIN = 1, function(ind) pmax(M[,ind[1]],M[,ind[2]],M[,ind[3]]))
  nrow(M) %in% colSums(comp_value)
  if(nrow(M) %in% colSums(comp_value)) {return(3)}  
  #Test for dominating number = 4
  combs <- t(combn(x = 1:ncol(M),m = 4))
  comp_value <- apply(X = combs, MARGIN = 1, function(ind) pmax(M[,ind[1]],M[,ind[2]],M[,ind[3]],M[,ind[4]]))
  nrow(M) %in% colSums(comp_value)
  if(nrow(M) %in% colSums(comp_value)) {return(4)}  
  #Test for dominating number = 5
  combs <- t(combn(x = 1:ncol(M),m = 5))
  comp_value <- apply(X = combs, MARGIN = 1, function(ind) pmax(M[,ind[1]],M[,ind[2]],M[,ind[3]],M[,ind[4]],M[,ind[5]]))
  nrow(M) %in% colSums(comp_value)
  if(nrow(M) %in% colSums(comp_value)) {return(5)} 
  #Test for dominating number = 6
  combs <- t(combn(x = 1:ncol(M),m = 6))
  comp_value <- apply(X = combs, MARGIN = 1, function(ind) pmax(M[,ind[1]],M[,ind[2]],M[,ind[3]],M[,ind[4]],M[,ind[5]],M[,ind[6]]))
  nrow(M) %in% colSums(comp_value)
  if(nrow(M) %in% colSums(comp_value)) {return(6)}   
  #Otherwise use a greedy set cover algorithm
  library(RcppGreedySetCover)
  S <- setsystem_from_matrix(M)
  invisible(capture.output(res <- greedySetCover(S,TRUE)))
  D <- uniqueN(res[,1])
  return(D)
}

maximal_points_from_matrix <- function(M,RowsOrCols = "Cols"){
  #Identify maximal points from a matrix
  if(RowsOrCols == "Cols") {
    z <- unlist(lapply(X = 1:ncol(M),function(k) all(unlist(lapply(X = 1:ncol(M),
                                                                   function(n) any(M[,k,drop = FALSE] > M[,n,drop = FALSE]) | 
                                                                     all(M[,k,drop = FALSE] == M[,n,drop = FALSE]))))))
  }
  if(RowsOrCols == "Rows") {
    z <- unlist(lapply(X = 1:nrow(M),function(k) all(unlist(lapply(X = 1:nrow(M),
                                                                   function(n) any(M[k,,drop = FALSE] > M[n,,drop = FALSE]) | 
                                                                     all(M[k,,drop = FALSE] == M[n,,drop = FALSE]))))))
  }  
  return(z)
}

minimal_points_from_matrix <- function(M,RowsOrCols = "Rows"){
  #Identify maximal points from a matrix
  if(RowsOrCols == "Cols") {
    z <- unlist(lapply(X = 1:ncol(M),function(k) all(unlist(lapply(X = 1:ncol(M),
                                                                   function(n) any(M[,k,drop = FALSE] < M[,n,drop = FALSE]) | 
                                                                     all(M[,k,drop = FALSE] == M[,n,drop = FALSE]))))))
    return(z)
  }
  if(RowsOrCols == "Rows") {
    z <- unlist(lapply(X = 1:nrow(M),function(k) all(unlist(lapply(X = 1:nrow(M),
                                                                   function(n) any(M[k,,drop = FALSE] < M[n,,drop = FALSE]) | 
                                                                     all(M[k,,drop = FALSE] == M[n,,drop = FALSE]))))))
    return(z)
  }  
}

remove_twins <- function(M){
  #Remove twin points from a matrix
  A <- M
  row_dups <- duplicated(A)
  col_dups <- duplicated(t(A))
  A <- A[!row_dups,!col_dups, drop = FALSE]
  return(A)
}

skeleton_bimorphic_form <- function(M){
  A <- M
  maximal_points <- maximal_points_from_matrix(A,RowsOrCols = "Cols")
  minimal_points <- minimal_points_from_matrix(A,RowsOrCols = "Rows")
  while(sum(maximal_points) < length(maximal_points) || sum(minimal_points) < length(minimal_points)) {
    A <- A[which(minimal_points),which(maximal_points),drop = FALSE]
    maximal_points <- maximal_points_from_matrix(A,RowsOrCols = "Cols")
    minimal_points <- minimal_points_from_matrix(A,RowsOrCols = "Rows")
  }
  A <- remove_twins(A)
  return(A)
}

canonical_form <- function(M) {
  i <- nrow(M)
  j <- ncol(M)
  I <- rowSums(M)
  J <- colSums(M)
  if(identical(I,rep(1,i)) & identical(J,rep(1,j))) {
    return(diag(i))
  }
  else if (identical(I,rep(j-1,i)) & identical(J,rep(i-1,j))) {
    return(1-diag(i))
  }
  else {
    return(M[order(rowSums(-M)),order(colSums(-M)),drop = FALSE])
  }
}

get_morphisms_from_AtoB <- function(A,B) {
  library(R.utils)
  #Explicitly Test for Morphism
  Amin <- 1:nrow(A)
  Bmin <- 1:nrow(B)
  Apl <- 1:ncol(A)
  Bpl <- 1:ncol(B)
  
  phimin <- permutations(n = length(Amin),r = length(Bmin), v = Amin, repeats.allowed = TRUE)
  phipl  <- permutations(n = length(Bpl),r = length(Apl), v = Bpl, repeats.allowed = TRUE)
  
  rowsandcols <- unlist(lapply(X = 1:nrow(phimin), function(k) lapply(X = 1:nrow(phipl), 
                                                                      function(n) list(phimin[k,],phipl[n,]))),recursive = FALSE)
  
  M1 <- lapply(X = 1:length(rowsandcols), function(k) A[rowsandcols[[k]][[1]],])
  M2 <- lapply(X = 1:length(rowsandcols), function(k) B[,rowsandcols[[k]][[2]]])
  morphism <- unlist(lapply(X = 1:length(rowsandcols), function(k) all(M2[[k]] >= M1[[k]])))
  
  z <- rowsandcols[which(morphism)]
  return(z)
}

matequal <- function(A, B) {
  if(is.matrix(A) && is.matrix(B) && all(dim(A) == dim(B))){return(all(A == B))}
  else{return(FALSE)}
}

morphism_from_AtoB <- function(A,B) {
  #Simply the relations to the Skeleton Bimorphic Forms
  A <- skeleton_bimorphic_form(A)
  B <- skeleton_bimorphic_form(B)
  #Identity
  if(matequal(A,B)){return(TRUE)}
  #Calculate dominating numbers and check for morphism existence conditions
  dA <- dominating_number(A)
  #A morphisms onto everything when d(A) DNE
  if(is.na(dA)){return(TRUE)}
  dB <- dominating_number(B)
  #No relation morphisms onto B when d(B) DNE
  if(is.na(dB)){return(FALSE)}
  #All relations morphism onto B when d(B) = 1
  if(dB == 1){return(TRUE)} 
  #d(A) < d(B) violates dominating number
  if(dA < dB){return(FALSE)} 
  #Duals
  Adual <- dual_relation(A)
  Bdual <- dual_relation(B)
  dA_dual <- dominating_number(Adual)
  dB_dual <- dominating_number(Bdual)
  #d(A_dual) > d(B_dual) violates dominating number
  if(dA_dual > dB_dual){return(FALSE)} 
  #The Two Ladder is a special case. All relations A such that dA = dA_dual = 2 are bimorphic with the two ladder.
  #If dA = dA_dual = 2, convert to two ladder
  if(dA ==2 & dA_dual ==2){
    A <- matrix(c(1,0,0,1),nrow = 2,ncol = 2,byrow = TRUE)    
  }
  if(dB ==2 & dB_dual ==2){
    B <- matrix(c(1,0,0,1),nrow = 2,ncol = 2,byrow = TRUE) 
  }
  #Create vectors representing the sets A-, B-, A+, B+
  Amin <- 1:nrow(A)
  Bmin <- 1:nrow(B)
  Apl <- 1:ncol(A)
  Bpl <- 1:ncol(B)
  
  #Look for obvious morphism existence conditions
  #Ladders morphisms onto ladders of lesser size
  if(dA == length(Apl) & length(Apl) == length(Amin) & dA > dB & dB == length(Bpl) & length(Bpl) == length(Bmin)){return(TRUE)} 
  #(Dual) Ladders morphisms onto ladders of lesser size
  if(dA_dual == length(Apl) & length(Apl) == length(Amin) & dA_dual < dB_dual & dB_dual == length(Bpl) & length(Bpl) == length(Bmin)){return(TRUE)} 
  #A ladder morphisms onto all relations with the same dominating number
  if(dA == length(Apl) & length(Apl) == length(Amin) & dA == dB){return(TRUE)} 
  #(Dual) A ladder morphisms onto all relations with the same dominating number
  if(dB_dual == length(Bpl) & length(Bpl) == length(Bmin) & dA_dual == dB_dual){return(TRUE)} 
  
  #If none of the above conditions apply, explicitly test for morphism
  #Define a recursive function for finding a viable phi+ given a phi-. This avoids the need to check all possible pairings. 
  phipl_recur <- function(phipl){
    #Check if a morphism is possible
    if(length(phipl) > 0){
      M1 <- A[phimin,1:length(phipl),drop = FALSE]
      M2 <- B[1:length(phimin),phipl,drop = FALSE]
      morphism_possible <- all(M2 >= M1)
      #If not possible, return FALSE
      if(!morphism_possible){
        return(FALSE)
      }
      #If is possible, and the function is completely defined, return TRUE
      if(length(phipl) == length(Apl) & morphism_possible){
        return(TRUE)
      }
    }
    #If it is possible, but the function isn't completely defined, continue recursively
    phipl_prev <- phipl
    phipl_name <- deparse(substitute(phipl))
    return_phrase <- ""
    #Test all possible "one-extensions" of the existing phi+
    for(i in Bpl){
      phipl_new <- c(phipl_prev,i) 
      assign(paste(phipl_name,i,sep = ""),phipl_new)
      return_phrase <- paste(return_phrase,"phipl_recur(",phipl_name,i,") | ",sep = "")
    }
    return_phrase <- substring(return_phrase,1,nchar(return_phrase)-3)
    return(eval(parse(text = return_phrase)))
  }
  #Generate all possible phi-
  phimin_list <- permutations(n = length(Amin),r = length(Bmin), v = Amin, repeats.allowed = TRUE)
  #Throw out any that don't have at least dA_dual distinct elements, since these are guaranteed to not result in morphism
  phimincount <- apply(X = phimin_list, MARGIN = 1, function(x) length(unique(x)))
  phimin_list <- phimin_list[which(phimincount >= dA_dual),]
  #Iterate over all possible phi-, use the phipl_recur function to look for viable phi+. 
  #If a morphism is found, stop searching and return TRUE
  for(i in 1:nrow(phimin_list)){
    phimin <- phimin_list[i,]
    root <- c()
    if(phipl_recur(root)){
      return(TRUE)
    }
  }
  #Otherwise, return FALSE
  return(FALSE)
} 

#Save Workspace
save.image(file = "TukeyMorphism_functionsLoaded.RData")

#####################################################################
################# End Function Definition ###########################
#####################################################################


#####################################################################
#################### Begin Relation Generation ######################
#####################################################################

Usize = 6 #number of elements in A-
Vsize = 6 #number of elements in A+

#This portion can be memory intensive, based on the sizes of A- and A+.
#Double check your capacity before running anything larger than 5x5.
#For example, 6x6 generates 109,453,344 distinct graphs, compared to 324,632 in the 5x5 case. 
#Use the function below to test sizes before running the rest of the code. 
count_of_graph_generation(Usize,Vsize)

#Due to memory constraints, generate in a loop and just keep the new skeleton bimorphic forms
loopcounter = 1
loopsize = 1000000 #Define loop size

#Generate finite relations of a given size
#Create all possible binary vectors of length j
#2^j possibilities at each step
perm <- permutations(n= 2, r = Vsize, v = 0:1, repeats.allowed = TRUE)

#Get "all" ixj matrices by combining (with replacement) the vectors from the previous step i times
#The zero vectors can be left out, since we are only interested in relations that have dominating families 
#((2^j-1)+i-1)!/(i!((2^j-1)-1)!) possibilities
combi <- combinations(n = nrow(perm)-1, r = Usize,v = 2:nrow(perm), repeats.allowed = TRUE)

#Initialize lists to store the skeleton matrices and graphs
#Include the empty relation
unique_skeleton_matrices <- list(matrix(0))
unique_skeleton_graphs <- list(graph_from_matrix(unique_skeleton_matrices[[1]],plotGraph = F))

#Begin Loop
while(loopcounter < nrow(combi)){
  a <- loopcounter
  b <- min(loopcounter+loopsize-1,nrow(combi))
  
  matrices <- lapply(X = a:b,function(n) matrix(perm[combi[n,],],byrow = F, nrow = Usize, ncol = Vsize))

  #Create Minimal Relations
  #Get the canonical form (order by degree) to make the duplicates more obvious
  skeleton_matrices <- lapply(X = matrices,function(M) canonical_form(skeleton_bimorphic_form(M)))
  
  #Find Unique (up to isomorphism) Minimal Relations
  #Remove exact duplicates
  unique_skeleton_matrices_stage <- unique(skeleton_matrices)
  unique_skeleton_graphs_stage <- lapply(X = unique_skeleton_matrices_stage,function(M) graph_from_matrix(M,plotGraph = F)) 
  
  unique_skeleton_matrices_stage <- append(unique_skeleton_matrices,unique_skeleton_matrices_stage)
  unique_skeleton_graphs_stage <- append(unique_skeleton_graphs,unique_skeleton_graphs_stage)
  
  #Find all graphs that are isomorphic
  isomorphic_list <- lapply(X = unique_skeleton_graphs_stage, 
                            function(G1) unlist(lapply(X = unique_skeleton_graphs_stage,
                                                       function(G2) isomorphic(G1,G2,method = 'vf2'))))
  isomorphism_classes <- unique(lapply(X = isomorphic_list, function(L) which(L)))
  
  #Find all matrices that are symmetric
  symmetric_list <- unlist(lapply(X = unique_skeleton_matrices_stage, function(M) t(isSymmetric(M))))
  symmetric_indices <- which(symmetric_list)
  
  #Choose a representative for each isomorphic class. Choose a symmetric relation if one exists.
  unique_skeleton_relations_final_indices <- unlist(lapply(X = isomorphism_classes, 
                                                          function(C){if(length(intersect(C,symmetric_indices)) == 0){C[1]}
                                                            else{intersect(C,symmetric_indices)[1]}}))
  
  unique_skeleton_matrices <- unique_skeleton_matrices_stage[unique_skeleton_relations_final_indices]
  unique_skeleton_graphs <- unique_skeleton_graphs_stage[unique_skeleton_relations_final_indices]
  
  save.image(file = "TukeyMorphism_inGraphGenerationLoop.RData")
  
  loopcounter <- b+1
} 
#End Loop

#####################################################################
#################### End Relation Generation ########################
#####################################################################


###################################################################################
### Create Data Frame for Storing Pairs of Relations and Checking for Morphisms ###
###################################################################################

#Check that minimal relations are closed under dual
unique_skeleton_matrices_duals <- lapply(X = unique_skeleton_matrices, function(M) dual_relation(M))
unique_skeleton_graphs_duals <- lapply(X = unique_skeleton_matrices_duals, function(M) graph_from_matrix(M,plotGraph = F)) 
isomorphic_list_duals <- lapply(X = unique_skeleton_graphs_duals, 
                                function(G1) unlist(lapply(X = unique_skeleton_graphs, 
                                                           function(G2) isomorphic(G1,G2,method = 'vf2'))))

#Calculate the Dominating Numbers, Dual Dominating Numbers, Size of A- and A+
Aminus <- unlist(lapply(X = unique_skeleton_matrices,function(M) nrow(M)))
Aplus  <- unlist(lapply(X = unique_skeleton_matrices,function(M) ncol(M)))
D      <- unlist(lapply(X = unique_skeleton_matrices,function(M) dominating_number(M)))
Ddual  <- unlist(lapply(X = unique_skeleton_matrices,function(M) dominating_number(dual_relation(M))))
Aindex <- 1:length(unique_skeleton_matrices)
dualIndex <- unlist(lapply(X = isomorphic_list_duals, function(L) if(length(which(L)) == 0){NA}else{which(L)}))
skeleton_characteristics <- cbind(Aindex,Aminus,Aplus,D,Ddual,dualIndex)
colnames(skeleton_characteristics) <- c("Aindex","Aminus","Aplus","d(A)","d(Adual)","dualIndex")

#Classification of morphisms for graphs
#Create a data frame with all combinations of minimal relations
morphisms <- merge(data.frame(skeleton_characteristics = skeleton_characteristics), 
                               data.frame(skeleton_characteristics = skeleton_characteristics), 
                               by = NULL)
names(morphisms) <- c("Aindex","Aminus","Aplus","dA","dA_dual","A_dualIndex",
                                  "Bindex","Bminus","Bplus","dB","dB_dual","B_dualIndex")


######################################################################
###Begin Classifying Morphism Existence for each Pair of Relations ###
######################################################################

#Classify based on existing lemmas
morphisms$MorphismFromAtoB <- with(morphisms, 
                                    ifelse(Aindex == Bindex,
                                            TRUE,  
                                    ifelse(is.na(dA),
                                            TRUE,       
                                    ifelse(is.na(dB),
                                            FALSE,
                                    ifelse(dB == 1,
                                            TRUE,
                                    ifelse(dA < dB,
                                            FALSE,
                                    ifelse(dA_dual > dB_dual,
                                            FALSE,   
                                    ifelse(dA == Aplus & Aplus == Aminus & dA > dB & dB == Bplus & Bplus == Bminus,
                                            TRUE,
                                    ifelse(dA_dual == Aplus & Aplus == Aminus & dA_dual < dB_dual & dB_dual == Bplus & Bplus == Bminus,
                                            TRUE,                                                         
                                    ifelse(dA == Aplus & Aplus == Aminus & dA == dB,
                                            TRUE,
                                    ifelse(dB_dual == Bplus & Bplus == Bminus & dA_dual == dB_dual,
                                            TRUE,
                                            ""
                                    )))))))))))

table1 <- table(morphisms$MorphismFromAtoB)
table1

#Add in "plus one" inclusions (i.e. A -> A' if one edge is added)
for(x in 1:length(unique_skeleton_matrices)) {
  X <- unique_skeleton_matrices[[x]]
  loopvector <- which(X == 0)
  E <- matrix(0,nrow = nrow(X),ncol = ncol(X))
  XPlusOne_matrices <- list()
  XPlusOne_graphs <- list()
  a <- 1
  for (i in loopvector) {
    EPlusOne <- E 
    EPlusOne[i] <- 1   
    XPlusOne_matrices[[a]] <- X + EPlusOne
    XPlusOne_graphs[[a]] <- graph_from_matrix(XPlusOne_matrices[[a]],plotGraph = F)
    a <- a+1
  }
  PlusOne_isomorphic_list <- lapply(X = XPlusOne_graphs, 
                            function(G1) unlist(lapply(X = unique_skeleton_graphs, function(G2) isomorphic(G1,G2,method = 'vf2'))))
  PlusOne_isomorphism_classes <- unlist(unique(lapply(X = PlusOne_isomorphic_list, function(L) which(L))))
  morphisms$MorphismFromAtoB[which(morphisms$Aindex == x & morphisms$Bindex %in% PlusOne_isomorphism_classes)] <- TRUE
}

table2 <- table(morphisms$MorphismFromAtoB)
table2

#Loop through Disjoint Unions, Dual Closure, and Transitive Closure as long as they are adding new morphisms
  #Disjoint Unions
  #Create all possible disjoint unions (that stay within the size limit)
  disjointUnions_stage <- list()
  for(i in 1:nrow(morphisms)){
    disjointUnions_stage[[i]]<- if(sum(morphisms[i,c("Aminus","Bminus")])<=Usize & sum(morphisms[i,c("Aplus","Bplus")])<=Vsize 
                                   & morphisms[i,"Aindex"] != 1 & morphisms[i,"Bindex"] != 1) {
      unique_skeleton_graphs[[morphisms[i,"Aindex"]]] %du% unique_skeleton_graphs[[morphisms[i,"Bindex"]]]}
    else{NULL}
  }
  disjointUnions <- disjointUnions_stage[unlist(lapply(disjointUnions_stage,function(k) !is.null(k)))]
  #Keep track of which two sub-relations created the union
  disjointUnions_indices<-morphisms[which(unlist(lapply(disjointUnions_stage,function(k) !is.null(k)))),c("Aindex","Bindex")]
  #See which relation each union is isomorphic to
  disjointUnions_isomorphic_list <- unlist(lapply(X = lapply(X = disjointUnions, function(G1) 
    unlist(lapply(X = unique_skeleton_graphs, function(G2) 
      isomorphic(G1,G2,method = 'vf2')))), function(L) which(L)))

#Counters
c <- 1
d <- 0

#Begin Loop
while (c > d) {
  #Limit to morphisms that involve these relations and haven't been decided
  morphisms_disjointUnions <- morphisms[which(morphisms$Aindex %in% disjointUnions_isomorphic_list & morphisms$MorphismFromAtoB == ''),]
  #Loop through and assign decision where possible based on the component pieces
  for(i in 1:nrow(morphisms_disjointUnions)){
    #Check if A1 -> B OR A2 -> B
    if(
       morphisms[which(morphisms$Aindex == disjointUnions_indices[min(which(disjointUnions_isomorphic_list == morphisms_disjointUnions[i,"Aindex"])),1] & morphisms$Bindex == morphisms_disjointUnions[i,"Bindex"]),"MorphismFromAtoB"] == "TRUE" | 
       morphisms[which(morphisms$Aindex == disjointUnions_indices[min(which(disjointUnions_isomorphic_list == morphisms_disjointUnions[i,"Aindex"])),2] & morphisms$Bindex == morphisms_disjointUnions[i,"Bindex"]),"MorphismFromAtoB"] == "TRUE"
    ){morphisms_disjointUnions[i,"MorphismFromAtoB"] <- "TRUE"}
    #Check if (A1 -> B1 AND A2 -> B2) OR (A1 -> B2 AND A2 -> B1)
    else if(morphisms_disjointUnions[i,"Bindex"] %in% disjointUnions_isomorphic_list){
      if (
          (
          morphisms[which(morphisms$Aindex == disjointUnions_indices[min(which(disjointUnions_isomorphic_list == morphisms_disjointUnions[i,"Aindex"])),1] & morphisms$Bindex == disjointUnions_indices[min(which(disjointUnions_isomorphic_list == morphisms_disjointUnions[i,"Bindex"])),1]),"MorphismFromAtoB"] == "TRUE" &
          morphisms[which(morphisms$Aindex == disjointUnions_indices[min(which(disjointUnions_isomorphic_list == morphisms_disjointUnions[i,"Aindex"])),2] & morphisms$Bindex == disjointUnions_indices[min(which(disjointUnions_isomorphic_list == morphisms_disjointUnions[i,"Bindex"])),2]),"MorphismFromAtoB"] == "TRUE"
          )
          |
          (
            morphisms[which(morphisms$Aindex == disjointUnions_indices[min(which(disjointUnions_isomorphic_list == morphisms_disjointUnions[i,"Aindex"])),1] & morphisms$Bindex == disjointUnions_indices[min(which(disjointUnions_isomorphic_list == morphisms_disjointUnions[i,"Bindex"])),2]),"MorphismFromAtoB"] == "TRUE" &
            morphisms[which(morphisms$Aindex == disjointUnions_indices[min(which(disjointUnions_isomorphic_list == morphisms_disjointUnions[i,"Aindex"])),2] & morphisms$Bindex == disjointUnions_indices[min(which(disjointUnions_isomorphic_list == morphisms_disjointUnions[i,"Bindex"])),1]),"MorphismFromAtoB"] == "TRUE"
          )
        )
        {morphisms_disjointUnions[i,"MorphismFromAtoB"] <- "TRUE"}
    }
  }
  #Add result back into table
  morphisms[which(morphisms$Aindex %in% disjointUnions_isomorphic_list & morphisms$MorphismFromAtoB == ''),"MorphismFromAtoB"] <- morphisms_disjointUnions[,"MorphismFromAtoB"]

  #Counter
  c <- table(morphisms$MorphismFromAtoB)[1]     
  
  #Dual Closure. If Bdual can morphism onto Adual, then A can morphism onto B, and vice versa
  morphism_sets_from_duals <- lapply(X = 1:length(unique_skeleton_matrices), 
                                     function(k) morphisms$A_dualIndex[which(morphisms$B_dualIndex == k & morphisms$MorphismFromAtoB == TRUE)])
  #Add dual closure back to the data frame
  for(k in 1:length(morphism_sets_from_duals)) {
    morphisms$MorphismFromAtoB[which(morphisms$Aindex == k & morphisms$Bindex %in% morphism_sets_from_duals[[k]])] <- TRUE
  }
  
  #Dual Closure. If Bdual can't morphism onto Adual, then A can't morphism onto B, and vice versa  
  non_morphism_sets_from_duals <- lapply(X = 1:length(unique_skeleton_matrices), 
                                         function(k) morphisms$A_dualIndex[which(morphisms$B_dualIndex == k & morphisms$MorphismFromAtoB == FALSE)])
  #Add dual closure back to the data frame
  for(k in 1:length(non_morphism_sets_from_duals)) {
      morphisms$MorphismFromAtoB[which(morphisms$Aindex == k & morphisms$Bindex %in% non_morphism_sets_from_duals[[k]])] <- FALSE
  }
  
  #Transitive Closure
  #List all B that can be morphismed onto by A
  morphism_sets <- lapply(X = 1:length(unique_skeleton_matrices),
                          function(k) morphisms$Bindex[which(morphisms$Aindex == k & morphisms$MorphismFromAtoB == TRUE)])
  #By transitivity, A can morphism onto anything that B can. Add all of those to the list for A.
  morphism_sets_with_transitivity <- lapply(X = 1:length(morphism_sets), function(k) 
    sort(unique(unlist(lapply(X = morphism_sets[[k]], function(n) morphism_sets[[n]])))))
  #Add transitive closures back to the data frame
  for(k in 1:length(morphism_sets_with_transitivity)) {
      morphisms$MorphismFromAtoB[which(morphisms$Aindex == k & morphisms$Bindex %in% morphism_sets_with_transitivity[[k]])] <- TRUE
  }
  
  #Contrapositive Transitive Closure
  #List all A that morphism onto B
  morphism_sets_onto <- lapply(X = 1:length(unique_skeleton_matrices),
                               function(k) morphisms$Aindex[which(morphisms$Bindex == k & morphisms$MorphismFromAtoB == TRUE)])
  non_morphism_sets <- lapply(X = 1:length(unique_skeleton_matrices),
                              function(k) morphisms$Bindex[which(morphisms$Aindex == k & morphisms$MorphismFromAtoB == FALSE)])
  #By the contrapositive of transitivity, If A -> B and not A->C, then not B->C.
  non_morphism_sets_with_transitivity <- lapply(X = 1:length(morphism_sets_onto), function(k) 
    sort(unique(unlist(lapply(X = morphism_sets_onto[[k]], function(n) non_morphism_sets[[n]])))))
  #Add transitive closures back to the data frame
  for(k in 1:length(non_morphism_sets_with_transitivity)) {
    morphisms$MorphismFromAtoB[which(morphisms$Aindex == k & morphisms$Bindex %in% non_morphism_sets_with_transitivity[[k]])] <- FALSE
  }
  
  #Counter
  d <- table(morphisms$MorphismFromAtoB)[1]
}
#End Loop

table3 <- table(morphisms$MorphismFromAtoB)
table3

#Looks for non-surjective phi_-
morphisms_nonsurjective <- morphisms[which(morphisms$MorphismFromAtoB == '' & morphisms$Aminus > morphisms$Bminus),]
for (j in 1:nrow(morphisms_nonsurjective)){
  A <- unique_skeleton_matrices[[morphisms_nonsurjective$Aindex[j]]]
  B <- unique_skeleton_matrices[[morphisms_nonsurjective$Bindex[j]]]
  N <- morphisms_nonsurjective$Aminus[j]
  R <- morphisms_nonsurjective$Bminus[j]
  C <- combinations(n = N,r = R)
  AprimeMorphism <- c()
  for (i in 1:nrow(C)) {
    Aprime <- A[C[i,],]
    AprimeMorphism[i] <- morphism_from_AtoB(Aprime,B)
  }
  #A -> B iff exists Aprime s.t. Aprime -> B
  morphisms_nonsurjective[j,"MorphismFromAtoB"] <- any(AprimeMorphism)
}
#Add result back into table
morphisms[which(morphisms$MorphismFromAtoB == '' & morphisms$Aminus > morphisms$Bminus),"MorphismFromAtoB"] <- morphisms_nonsurjective[,"MorphismFromAtoB"]

#Counters
e <- 1
f <- 0

#Begin Loop
while (e > f) {
  #Dual Closure. If Bdual can morphism onto Adual, then A can morphism onto B, and vice versa
  morphism_sets_from_duals <- lapply(X = 1:length(unique_skeleton_matrices), 
                                     function(k) morphisms$A_dualIndex[which(morphisms$B_dualIndex == k & morphisms$MorphismFromAtoB == TRUE)])
  #Add dual closure back to the data frame
  for(k in 1:length(morphism_sets_from_duals)) {
    morphisms$MorphismFromAtoB[which(morphisms$Aindex == k & morphisms$Bindex %in% morphism_sets_from_duals[[k]])] <- TRUE
  }
  
  #Dual Closure. If Bdual can't morphism onto Adual, then A can't morphism onto B, and vice versa  
  non_morphism_sets_from_duals <- lapply(X = 1:length(unique_skeleton_matrices), 
                                         function(k) morphisms$A_dualIndex[which(morphisms$B_dualIndex == k & morphisms$MorphismFromAtoB == FALSE)])
  #Add dual closure back to the data frame
  for(k in 1:length(non_morphism_sets_from_duals)) {
    morphisms$MorphismFromAtoB[which(morphisms$Aindex == k & morphisms$Bindex %in% non_morphism_sets_from_duals[[k]])] <- FALSE
  }
  
  #Counter
  e <- table(morphisms$MorphismFromAtoB)[1]  
  
  #Transitive Closure
  #List all B that can be morphismed onto by A
  morphism_sets <- lapply(X = 1:length(unique_skeleton_matrices),
                          function(k) morphisms$Bindex[which(morphisms$Aindex == k & morphisms$MorphismFromAtoB == TRUE)])
  #By transitivity, A can morphism onto anything that B can. Add all of those to the list for A.
  morphism_sets_with_transitivity <- lapply(X = 1:length(morphism_sets), function(k) 
    sort(unique(unlist(lapply(X = morphism_sets[[k]], function(n) morphism_sets[[n]])))))
  #Add transitive closures back to the data frame
  for(k in 1:length(morphism_sets_with_transitivity)) {
    morphisms$MorphismFromAtoB[which(morphisms$Aindex == k & morphisms$Bindex %in% morphism_sets_with_transitivity[[k]])] <- TRUE
  }
  
  #Contrapositive Transitive Closure
  #List all A that morphism onto B
  morphism_sets_onto <- lapply(X = 1:length(unique_skeleton_matrices),
                               function(k) morphisms$Aindex[which(morphisms$Bindex == k & morphisms$MorphismFromAtoB == TRUE)])
  non_morphism_sets <- lapply(X = 1:length(unique_skeleton_matrices),
                              function(k) morphisms$Bindex[which(morphisms$Aindex == k & morphisms$MorphismFromAtoB == FALSE)])
  #By the contrapositive of transitivity, If A -> B and not A->C, then not B->C.
  non_morphism_sets_with_transitivity <- lapply(X = 1:length(morphism_sets_onto), function(k) 
    sort(unique(unlist(lapply(X = morphism_sets_onto[[k]], function(n) non_morphism_sets[[n]])))))
  #Add transitive closures back to the data frame
  for(k in 1:length(non_morphism_sets_with_transitivity)) {
    morphisms$MorphismFromAtoB[which(morphisms$Aindex == k & morphisms$Bindex %in% non_morphism_sets_with_transitivity[[k]])] <- FALSE
  }
  
  #Counter
  f <- table(morphisms$MorphismFromAtoB)[1]
}
#End Loop

table4 <- table(morphisms$MorphismFromAtoB)
table4

#See what % of pairs have been classified so far
#For the 6x6 case, it is 98.6%
(1-table4[1]/sum(table4))*100

#Save Workspace
save.image(file = "TukeyMorphism_readyToCheckMorphisms.RData")

#Explicitly check for morphisms for remaining pairs
primary_morphisms <- sqldf("select a.* 
                           from morphisms a 
                           join morphisms b 
                            on a.B_dualIndex = b.Aindex 
                            and a.A_dualIndex = b.Bindex
                           where a.Aindex < b.Aindex --limit to one pair A->B
                           and a.MorphismFromAtoB = '' --morphism hasn't been determined")

MorphismFromAtoB <- unlist(lapply(X = 1:nrow(primary_morphisms)
                            , function(k) morphism_from_AtoB(unique_skeleton_matrices[[primary_morphisms[k,1]]],
                                                             unique_skeleton_matrices[[primary_morphisms[k,7]]])))

primary_morphisms <- cbind(primary_morphisms[,-13],MorphismFromAtoB)

#Add explicit morphisms into data
morphisms_final <- merge(morphisms, primary_morphisms, by = c("Aindex","Bindex"), all.x = TRUE, suffixes = c("",".pm"))
morphisms_final <- morphisms_final[,c(names(morphisms),"MorphismFromAtoB.pm")]
morphisms_final$MorphismFromAtoB_final <- with(morphisms_final, ifelse(is.na(MorphismFromAtoB.pm),MorphismFromAtoB,MorphismFromAtoB.pm))

#Dual Closure. If Bdual can morphism onto Adual, then A can morphism onto B, and vice versa
morphism_sets_from_duals <- lapply(X = 1:length(unique_skeleton_matrices), 
                                   function(k) morphisms_final$A_dualIndex[which(morphisms_final$B_dualIndex == k & morphisms_final$MorphismFromAtoB_final == TRUE)])
#Add dual closure back to the data frame
for(k in 1:length(morphism_sets_from_duals)) {
  morphisms_final$MorphismFromAtoB_final[which(morphisms_final$Aindex == k & morphisms_final$Bindex %in% morphism_sets_from_duals[[k]])] <- TRUE
}

#Dual Closure. If Bdual can't morphism onto Adual, then A can't morphism onto B, and vice versa  
non_morphism_sets_from_duals <- lapply(X = 1:length(unique_skeleton_matrices), 
                                       function(k) morphisms_final$A_dualIndex[which(morphisms_final$B_dualIndex == k & morphisms_final$MorphismFromAtoB_final == FALSE)])
#Add dual closure back to the data frame
for(k in 1:length(non_morphism_sets_from_duals)) {
  morphisms_final$MorphismFromAtoB_final[which(morphisms_final$Aindex == k & morphisms_final$Bindex %in% non_morphism_sets_from_duals[[k]])] <- FALSE
}

#See summary of final classification
table(morphisms_final$MorphismFromAtoB_final)


#Save Workspace
save.image(file = "TukeyMorphism_final.RData")


#################################################
##################### Output ####################
#################################################

#Convert to Logical Matrix
morphism_matrix_logical <- matrix('',nrow = nrow(skeleton_characteristics),ncol = nrow(skeleton_characteristics))

for(i in 1:nrow(skeleton_characteristics)){
  for(j in 1:nrow(skeleton_characteristics)){
    morphism_matrix_logical[i,j] <- if(morphisms_final[which(morphisms_final$Aindex == i & morphisms_final$Bindex == j),'MorphismFromAtoB_final'] == "TRUE"){TRUE}
    else if(morphisms_final[which(morphisms_final$Aindex == i & morphisms_final$Bindex == j),'MorphismFromAtoB_final'] == "FALSE"){FALSE}
  }
}

morphism_matrix_logical <- matrix(as.logical(morphism_matrix_logical),ncol = ncol(morphism_matrix_logical))

#Print out Morphism Matrix
write.csv(morphism_matrix_logical,"Morphisms_6x6.csv") 

#Print out Morphism Characteristics 
write.csv(skeleton_characteristics,"Skeleton Characteristics_6x6.csv")

#Print out images of the graphs
#library(grDevices)
#lapply(X = 1:length(unique_skeleton_matrices),
#       function(k) {
#         M <- unique_skeleton_matrices[[k]]
#         mypath <- file.path("GraphicalOutput/",paste0("[",k,"].jpg", sep = ""))
#         jpeg(file=mypath)
#         graph_from_matrix(M, plotGraph = T,
#                           graphName = "",
#                           saveGraph = F)
#         dev.off()
#       })

#Hasse Diagram
#hasse(morphism_matrix_logical,parameters = list(clusterMerge = TRUE))

#Save Workspace
save.image(file = "TukeyMorphism_final.RData")
