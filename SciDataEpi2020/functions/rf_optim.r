
#   RFOPTIM ~ ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##   
			rf.optim <- function(counts.matrix, classes, iterations){
			iterations = iterations
			pb <- txtProgressBar(max=iterations, style=3)
			rf.optim.res <- list()
			for(i in 1:iterations){
				sample.vec <- partition.samples(y, 0.6)
#				train.x <- counts.matrix[sample.vec==1]
				train.y <- classes[sample.vec==1]
#				test.x <- counts.matrix[sample.vec==2]
				test.y <- classes[sample.vec==2]


				CM2train <- counts.matrix[sample.vec==1,]
				CM2test <-  counts.matrix[sample.vec==2,]


				rf <- randomForest(x=CM2train, y=train.y, ntree=1000)
				RF1.error <- confusionMatrix(data=predict(rf,CM2train), ref=train.y)
				RF2.error <- confusionMatrix(data=predict(rf,CM2test), ref=test.y)
			

				#rf <- randomForest(x=counts.matrix[sample.vec==1,feature.vector], y=train.y, ntree=100)
			
			#RF1.error <- confusionMatrix(data=predict(rf,counts.matrix[sample.vec==1,feature.vector]), ref=train.y)
			#RF2.error <- confusionMatrix(data=predict(rf,counts.matrix[sample.vec==2,feature.vector]), ref=test.y)
			setTxtProgressBar(pb, i)
			rf.optim.res[[i]] <- list("rf"= rf, "RFtrain.error" = RF1.error, "RFtest.error" = RF2.error)
			
			}
			return(rf.optim.res)
			}

#	rf.optim.res <- rf.optim(x,y,1000)
#	
#	n.genes2keep = 100
#  important.genes <- data.frame(table(unlist(lapply(1:length(rf.optim.res), function(i) {order(-rf.optim.res[[i]]$rf$importance[,1])[1:n.genes2keep]}))))
#  important.genes.vector <- as.integer(as.character(important.genes[,1]))
#  re.order <- order(-important.genes[,2])
#	important.genes <- important.genes[re.order,]
#	important.genes[,1] <- rownames(rf.optim.res[[1]]$rf$importance)[as.integer(as.character(important.genes[,1]))]

#	importance.scores.mean <- sapply(1:length(rf.optim.res), function(i) {rf.optim.res[[i]]$rf$importance})
#	importance.scores.sd <- sapply(1:length(importance.scores.mean[,1]), function(i) {sd(importance.scores.mean[i,])})
#	importance.scores.mean <- rowMeans(importance.scores.mean)
#  importance.scores.mean <- importance.scores.mean[important.genes.vector[re.order]]

#	rf.df <- data.frame(as.character(features), importance.scores.mean, stringsAsFactors=F)
#	rf.df <- rf.df[order(-rf.df[,2]),]

#	evidence.rf <- rf.df[,1]

