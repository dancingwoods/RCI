# Bronwyn Woods                                               
# 2013                                                         
#                                                              
# Summary: GUI for segmentation.                                                
# 
# License information
# This file is part of the R package RCI.
# 
# RCI is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# RCI is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with RCI.  If not, see <http://www.gnu.org/licenses/>.
# 
###############################################################################

# Uses RefClass to return the values from the GUI.  The initializer fires up the GUI, and then
# the result is put into the maskmat field which can then be retrieved with
# foo = CellID$new(db)
# result = foo$maskmat
CellID <- suppressWarnings(setRefClass(
	"CellID",
	fields=list(
		maskmat = "ANY"
	),
	methods = list(
		initialize = function(db){
			# A matrix with a row for each mask with the mask id and the current label and current segmentation
			selmat <- dbGetQuery(db$db, "select id, label, segmentation from masks order by id")

			mimg1 <- GetImage(db, 'mimg1')
			mimg1eq <- GetImage(db, 'mimg1eq')
			mimg2 <- GetImage(db, 'mimg2')
			mimg2eq <- GetImage(db, 'mimg2eq')

			dims <- as.integer(dbGetQuery(db$db, "select nx, ny from experiment")[1,])
			nx <- dims[1]
			ny <- dims[2]
	
			maskmatrixA <- matrix(NA, ny, nx)
			maskmatrixB <- matrix(NA, ny, nx)
	
			ImgHandlerA <- function(h, ...){
				xval <- round(h$x*nx)[1]
				yval <- round((1-h$y)*ny)[1]
		
				if(is.na(maskmatrixA[yval, xval])){
					maskmatrixA[yval, xval] <<- 1
				}else{
					maskmatrixA[yval, xval] <<- NA
				}
	
				RedrawImages()
				SaveMasks()
			}
	
			ImgHandlerB <- function(h, ...){
				xval <- round(h$x*nx)[1]
				yval <- round((1-h$y)*ny)[1]
		
				if(is.na(maskmatrixA[yval, xval])){
					maskmatrixB[yval, xval] <<- 1
				}else{
					maskmatrixB[yval, xval] <<- NA
				}
	
				RedrawImages()
				SaveMasks()
			}
	
			SaveMasks <- function(){
				maskmat <<- matrix(0, ny, nx)
				maskmat[which(!is.na(maskmatrixA))] <<- 1
				maskmat[which(!is.na(maskmatrixB))] <<- 2
			}
	
			RedrawImages <- function(){
				visible(imgwindow1) <- TRUE			
				Image(mimg1, useRaster=T)
				PlotMask(maskmatrixA, rgb=c(0,0,1))
				PlotMask(maskmatrixB, rgb=c(0,1,0))
				visible(imgwindow2) <- TRUE			
				Image(mimg2, useRaster=T)
				PlotMask(maskmatrixA, rgb=c(0,0,1))
				PlotMask(maskmatrixB, rgb=c(0,1,0))
				visible(imgwindow3) <- TRUE			
				Image(mimg1eq, useRaster=T)
				PlotMask(maskmatrixA, rgb=c(0,0,1))
				PlotMask(maskmatrixB, rgb=c(0,1,0))
				visible(imgwindow4) <- TRUE						
				Image(mimg2eq, useRaster=T)
				PlotMask(maskmatrixA, rgb=c(0,0,1))
				PlotMask(maskmatrixB, rgb=c(0,1,0))
		
			}

			########## Layout
	
			##The widgets
			win <- gwindow("Calcium Imaging Segmenter", spacing=10)
			maingroup <- ggroup(horizontal=FALSE, cont=win, spacing=5)
			vexpgroup <- ggroup(horizontal=FALSE, cont=maingroup, expand=TRUE, spacing=10)
			controlgroup <- ggroup(horizontal=TRUE, cont=maingroup, expand=TRUE, spacing=10)	

			############################
			##### view experiment group
		
			imgsubgroupA <- ggroup(cont=vexpgroup)
			imgsubgroupB <- ggroup(cont=vexpgroup)


			imgwindow1 <- ggraphics(container=imgsubgroupA)
			addHandlerChanged(imgwindow1, ImgHandlerA)
	
			imgwindow2 <- ggraphics(container=imgsubgroupA)
			addHandlerChanged(imgwindow2, ImgHandlerB)

			imgwindow3 <- ggraphics(container=imgsubgroupB)
			addHandlerChanged(imgwindow3, ImgHandlerA)

			imgwindow4 <- ggraphics(container=imgsubgroupB)
			addHandlerChanged(imgwindow4, ImgHandlerB)
	
			visible(imgwindow1) <- TRUE
			par(mar=c(0,0,0,0))
			Image(mimg1, useRaster=T)

			visible(imgwindow2) <- TRUE
			par(mar=c(0,0,0,0))
			Image(mimg2, useRaster=T)

			visible(imgwindow3) <- TRUE
			par(mar=c(0,0,0,0))
			Image(mimg1eq, useRaster=T)

			visible(imgwindow4) <- TRUE
			par(mar=c(0,0,0,0))
			Image(mimg2eq, useRaster=T)

		}
	)
))


# NOTE: this is currently broken in R 3.0 due to problems with RGtk2 - use old version of R
# Requires that mean images have been added to the database
#-
#' Opens the GUI viewer to manipulate the segmentation process.
#' 
#' @param db if specified, the viewer opens with the given dbController (looking
#' in the directories stored in that object)
#' @param cf a classifier object.  Needed to allow redoing segmentation after correcting labels
#' 
#' @return NULL
#' @export
#-
ViewCI <- function(db, cf=NULL){
	#options(guiToolkit="tctkl")
	
	# A matrix with a row for each mask with the mask id and the current label and current segmentation
	selmat <- dbGetQuery(db$db, "select id, label, segmentation from masks order by id")

	mimg1 <- GetImage(db, 'mimg1')
	mimg1eq <- GetImage(db, 'mimg1eq')
	mimg2 <- GetImage(db, 'mimg2')
	mimg2eq <- GetImage(db, 'mimg2eq')

	dims <- as.integer(dbGetQuery(db$db, "select nx, ny from experiment")[1,])
	nx <- dims[1]
	ny <- dims[2]
	
	# Just to keep track of the IDs
	negID <- 1
	neuronID <- 2
	astroID <- 3
	unknownID <- 0
	
	negColor <- c(1,0,0)
	neuronColor <- c(0,1,0)	
	astroColor <- c(0,0,1)
	unknownColor <- c(0.9,0.5,0.5)
	
	sms <- GetSparseMasks(db)

	# If pred=T, show predicted masks
	AddLabeledMasks <- function(pred=F){
		if(pred){
			wmasks <- 3
		}else{
			wmasks <- 2
		}
		
		visible(imgwindow1) <- TRUE			
		Image(mimg1, useRaster=T)
		visible(imgwindow2) <- TRUE			
		Image(mimg2, useRaster=T)
		visible(imgwindow3) <- TRUE			
		Image(mimg1eq, useRaster=T)
		visible(imgwindow4) <- TRUE						
		Image(mimg2eq, useRaster=T)
		
		# Plot neurons in green
		neuronmat = matrix(0, ny, nx)		
		for(i in 1:nrow(selmat)){
			if(!is.na(selmat[i,wmasks]) && abs(selmat[i,wmasks])==neuronID){
				neuronmat <- neuronmat+SparseToMatrix(sms[[i]], ny, nx, background=0)
			}
		}
		neuronmat[which(neuronmat>0)]=1
		if(max(neuronmat)>0){		
			visible(imgwindow1) <- TRUE
			PlotMaskSet(neuronmat, rgb=neuronColor, useRaster=T)
			visible(imgwindow2) <- TRUE
			PlotMaskSet(neuronmat, rgb=neuronColor, useRaster=T)
			visible(imgwindow3) <- TRUE
			PlotMaskSet(neuronmat, rgb=neuronColor, useRaster=T)
			visible(imgwindow4) <- TRUE
			PlotMaskSet(neuronmat, rgb=neuronColor, useRaster=T)
		}
		
		# Plot astrocytes in blue
		astromat = matrix(0, ny, nx)		
		for(i in 1:nrow(selmat)){
			if(!is.na(selmat[i,wmasks]) && abs(selmat[i,wmasks])==astroID){
				astromat <- astromat+SparseToMatrix(sms[[i]], ny, nx, background=0)
			}
		}
		astromat[which(astromat>0)]=1
		if(max(astromat)>0){		
			visible(imgwindow1) <- TRUE
			PlotMaskSet(astromat, rgb=astroColor, useRaster=T)
			visible(imgwindow2) <- TRUE
			PlotMaskSet(astromat, rgb=astroColor, useRaster=T)
			visible(imgwindow3) <- TRUE
			PlotMaskSet(astromat, rgb=astroColor, useRaster=T)
			visible(imgwindow4) <- TRUE
			PlotMaskSet(astromat, rgb=astroColor, useRaster=T)
		}
		
		
	}
	
	ToggleMasks <- function(h, ...){
		val <- svalue(togglemasks)
		if(val=="No Labels"){
			visible(imgwindow1) <- TRUE			
			Image(mimg1, useRaster=T)
			visible(imgwindow2) <- TRUE			
			Image(mimg2, useRaster=T)
			visible(imgwindow3) <- TRUE			
			Image(mimg1eq, useRaster=T)
			visible(imgwindow4) <- TRUE						
			Image(mimg2eq, useRaster=T)
		}else if(val=="Hand Labels"){
			AddLabeledMasks(F)
		}else{
			AddLabeledMasks(T)
		}
	}
			
	RedoSeg <- function(h, ...){
		print("Merging data")
		data <- cf$training
		tmp <- PullData(db)
		if(nrow(tmp)>0){
			if( ncol(tmp) == (ncol(data)-1) ){
				namecol <- as.data.frame(rep(dbGetQuery(db$db, "select tag from experiment")[1,1], nrow(tmp)))
				names(namecol) <- c("experiment")
				tmp <- cbind(namecol, tmp)
		 		data <- rbind(data, tmp)
			}
		}
		print("Retraining Classifier")
		cf <<- InitiateMaskClassifier(data)
		print("Recomputing Predictions")
		PredictExperiment(cf, db)
		selmat <<- dbGetQuery(db$db, "select id, label, segmentation from masks order by id")
		ToggleMasks(h)
	}

	InitiateFromClass <- function(h, ...){
		data <- PullData(db, F)
		data[is.na(data)]=0
		labs <- data[,which(names(data)=="label")]
		data  <- data[, -(which(names(data)=="label"))]
		classes <- predict(cf$cf, data)
		classes[which(classes==1)]=0
		
		res <- cbind(data[, which(names(data)=="id")], labs, classes)
		res <- res[order(res[,1]),]

		

		w <- c(which(is.na(res[,2])), which(res[,2]==0))
		res[w,2] <- res[w,3]
		
		res <- as.data.frame(res[,c(2,1)])
		names(res) <- c("label", "id")
		print(table(res[,1]))

		dbBeginTransaction(db$db)
		dbGetPreparedQuery(db$db, "UPDATE masks SET label=? where id=?", bind.data=res)
		dbCommit(db$db)
		selmat <<- dbGetQuery(db$db, "select id, label, segmentation from masks order by id")
	}

	# When someone clicks on or selects a region in the image
	ImgHandler <- function(h, ...){

		PlotPage <- function(pnum){
			for(i in 1:5){
				snum = 5*(pnum-1)+(i-1)
				for(j in 1:5){
					if(snum*5+j<=length(containedmasks)){
						id = -1* containedmasks[[sizeorder[snum*5+j]]][1]
						label = selmat[which(selmat[,1]==id),2]
						
						if(is.na(label)){
							rgb <- unknownColor
						}else if(label==negID){
							rgb <- negColor
						}else if(label==neuronID){
							rgb <- neuronColor
						}else if(label==astroID){
							rgb <- astroColor
						}else{
							rgb <- unknownColor
						}
						
						visible(ref.mswin[[(i-1)*5+j]]) <- TRUE
						par(mar=c(0.1,0.2,0.1,0.1), mfrow=c(1,2))
						Image(mimg1[yrange[1]:yrange[2], xrange[1]:xrange[2]], useRaster=T)
						PlotMask(SparseToMatrix(containedmasks[[sizeorder[snum*5+j]]], ny, nx)[yrange[1]:yrange[2], xrange[1]:xrange[2]], rgb=rgb, useRaster=T)
						par(mar=c(0.1,0.1,0.1,0.2))
						Image(mimg2[yrange[1]:yrange[2], xrange[1]:xrange[2]], useRaster=T)
						PlotMask(SparseToMatrix(containedmasks[[sizeorder[snum*5+j]]], ny, nx)[yrange[1]:yrange[2], xrange[1]:xrange[2]], rgb=rgb, useRaster=T)
					}else{
						visible(ref.mswin[[(i-1)*5+j]]) <- TRUE
						mat = matrix(1, diff(yrange), diff(xrange))
						mat[,1]=0
						Image(mat)
					}
				}
			}
			if(25*(pnum) >=length(containedmasks)){
				enabled(mswin.nextbutton) <- FALSE
			}else{
				enabled(mswin.nextbutton) <- TRUE
			}
			if(pnum>1){
				enabled(mswin.prevbutton) <- TRUE
			}else{
				enabled(mswin.prevbutton) <- FALSE
			}
		}			
		
		PagePrev <- function(h, ...){
			curpage <<- curpage-1
			PlotPage(curpage)
			svalue(mswin.pagelabel) <- paste("Showing ", (curpage-1)*25+1,"-", min(curpage*25,length(containedmasks)), " of ", length(containedmasks), sep="")				
		}
		
		PageNext <- function(h, ...){
			curpage <<- curpage+1
			PlotPage(curpage)
			svalue(mswin.pagelabel) <- paste("Showing ", (curpage-1)*25+1,"-", min(curpage*25,length(containedmasks)), " of ", length(containedmasks), sep="")
		}
		
		MarkNegative <- function(h, ...){
			# mark any unmarked masks in the current set as false positives
			for(i in 1:length(containedmasks)){
				selmatrow = which(selmat[,1]==-1*containedmasks[[i]][1])
				curlabel = selmat[selmatrow,2]
				if(is.na(curlabel) || curlabel==0){
					selmat[selmatrow,2] <<- negID
					SetMaskLabel(db, -1*containedmasks[[i]][1], negID)
				}
			}
			PlotPage(curpage)
		}
	
		
		CloseLabelWindow <- function(h, ...){
			if(svalue(togglemasks)=="Final Segmentation"){
				AddLabeledMasks(T)
			}else{
				AddLabeledMasks(F)
				svalue(togglemasks) <<- "Hand Labels"
			}
			dispose(mswin)
		}
	
		SetLabel <- function(index){
			ind = (curpage-1)*25+index
			if(ind<=length(containedmasks)){
				id = -1* containedmasks[[sizeorder[ind]]][1]
				label = selmat[which(selmat[,1]==id),2]
				if(is.na(label)){
					# Currently unlabeled
					# Change to neuron
					label <- neuronID
					rgb=neuronColor					
				}else if(label==negID){
					# Currently marked as negative
					# Change to unlabeled
					label <- unknownID
					rgb=unknownColor
				}else if(label==neuronID){
					# Currently marked as cell
					# change to astrocyte
					label <- astroID
					rgb=astroColor
				}else if(label==astroID){
					# Currently marked as astrocyte
					# change to negative
					label <- negID
					rgb <- negColor
				}else{
					# Currently unlabeled
					# Change to neuron
					label <- neuronID
					rgb=neuronColor
				}
				# change selmat
				selmat[which(selmat[,1]==id),2] <<- label
				# update in database
				SetMaskLabel(db, id, label)	
				
				visible(ref.mswin[[index]]) <- TRUE
				par(mar=c(0.1,0.2,0.1,0.1), mfrow=c(1,2))
				Image(mimg1[yrange[1]:yrange[2], xrange[1]:xrange[2]], useRaster=T)
				PlotMask(SparseToMatrix(containedmasks[[sizeorder[ind]]], ny, nx)[yrange[1]:yrange[2], xrange[1]:xrange[2]], rgb=rgb, useRaster=T)
				par(mar=c(0.1,0.1,0.1,0.2))
				Image(mimg2[yrange[1]:yrange[2], xrange[1]:xrange[2]], useRaster=T)
				PlotMask(SparseToMatrix(containedmasks[[sizeorder[ind]]], ny, nx)[yrange[1]:yrange[2], xrange[1]:xrange[2]], rgb=rgb, useRaster=T)
			}
		}
		
		SetLabel1 <- function(h, ...){SetLabel(1)}
		SetLabel2 <- function(h, ...){SetLabel(2)}
		SetLabel3 <- function(h, ...){SetLabel(3)}
		SetLabel4 <- function(h, ...){SetLabel(4)}
		SetLabel5 <- function(h, ...){SetLabel(5)}
		SetLabel6 <- function(h, ...){SetLabel(6)}
		SetLabel7 <- function(h, ...){SetLabel(7)}
		SetLabel8 <- function(h, ...){SetLabel(8)}
		SetLabel9 <- function(h, ...){SetLabel(9)}
		SetLabel10 <- function(h, ...){SetLabel(10)}
		SetLabel11 <- function(h, ...){SetLabel(11)}
		SetLabel12 <- function(h, ...){SetLabel(12)}
		SetLabel13 <- function(h, ...){SetLabel(13)}
		SetLabel14 <- function(h, ...){SetLabel(14)}
		SetLabel15 <- function(h, ...){SetLabel(15)}
		SetLabel16 <- function(h, ...){SetLabel(16)}
		SetLabel17 <- function(h, ...){SetLabel(17)}
		SetLabel18 <- function(h, ...){SetLabel(18)}
		SetLabel19 <- function(h, ...){SetLabel(19)}
		SetLabel20 <- function(h, ...){SetLabel(20)}
		SetLabel21 <- function(h, ...){SetLabel(21)}
		SetLabel22 <- function(h, ...){SetLabel(22)}
		SetLabel23 <- function(h, ...){SetLabel(23)}
		SetLabel24 <- function(h, ...){SetLabel(24)}
		SetLabel25 <- function(h, ...){SetLabel(25)}
		
		xrange <- round(h$x*nx)
		yrange <- round((1-h$y)*ny)[2:1]
	
		if(diff(xrange)==0 && diff(yrange)==0){
			#only one point selected
		}else{
			#rectangle selected
			region <- matrix(NA, ny, nx)
			region[yrange[1]:yrange[2], xrange[1]:xrange[2]] <- 1

			containedmasks <- GetInnerMasks(region, sms)
			sizeorder <- order(sapply(containedmasks, length), decreasing=T)
	
			# Create the window to allow users to label the masks in ovlp
			mswin <- gwindow("Label Candidate Masks", visible=FALSE)
			mswin.imggroup <- ggroup(container=mswin, horizontal=FALSE)
			mswin.refgroup <- ggroup(container=mswin.imggroup, expand=FALSE)
			mswin.junk1 <- ggraphics(width=400, height=100, container=mswin.refgroup)			
			mswin.refimg1 <- ggraphics(width=100, height=100, container=mswin.refgroup)
			mswin.refimg2 <- ggraphics(width=100, height=100, container=mswin.refgroup)
			mswin.junk2 <- ggraphics(width=400, height=100, container=mswin.refgroup)			
	
			mswin.imgrow1 <- ggroup(container=mswin.imggroup)
			mswin.imgrow2 <- ggroup(container=mswin.imggroup)
			mswin.imgrow3 <- ggroup(container=mswin.imggroup)
			mswin.imgrow4 <- ggroup(container=mswin.imggroup)
			mswin.imgrow5 <- ggroup(container=mswin.imggroup)
			
			imgmask1 <- ggraphics(width=200, height= 100, container=mswin.imgrow1)
			addHandlerChanged(imgmask1, SetLabel1)
			imgmask2 <- ggraphics(width=200, height= 100, container=mswin.imgrow1)
			addHandlerChanged(imgmask2, SetLabel2)
			imgmask3 <- ggraphics(width=200, height= 100, container=mswin.imgrow1)
			addHandlerChanged(imgmask3, SetLabel3)
			imgmask4 <- ggraphics(width=200, height= 100, container=mswin.imgrow1)
			addHandlerChanged(imgmask4, SetLabel4)
			imgmask5 <- ggraphics(width=200, height= 100, container=mswin.imgrow1)
			addHandlerChanged(imgmask5, SetLabel5)
			imgmask6 <- ggraphics(width=200, height= 100, container=mswin.imgrow2)
			addHandlerChanged(imgmask6, SetLabel6)
			imgmask7 <- ggraphics(width=200, height= 100, container=mswin.imgrow2)
			addHandlerChanged(imgmask7, SetLabel7)
			imgmask8 <- ggraphics(width=200, height= 100, container=mswin.imgrow2)
			addHandlerChanged(imgmask8, SetLabel8)
			imgmask9 <- ggraphics(width=200, height= 100, container=mswin.imgrow2)
			addHandlerChanged(imgmask9, SetLabel9)
			imgmask10 <- ggraphics(width=200, height= 100, container=mswin.imgrow2)
			addHandlerChanged(imgmask10, SetLabel10)
			imgmask11 <- ggraphics(width=200, height= 100, container=mswin.imgrow3)
			addHandlerChanged(imgmask11, SetLabel11)
			imgmask12 <- ggraphics(width=200, height= 100, container=mswin.imgrow3)
			addHandlerChanged(imgmask12, SetLabel12)
			imgmask13 <- ggraphics(width=200, height= 100, container=mswin.imgrow3)
			addHandlerChanged(imgmask13, SetLabel13)
			imgmask14 <- ggraphics(width=200, height= 100, container=mswin.imgrow3)
			addHandlerChanged(imgmask14, SetLabel14)
			imgmask15 <- ggraphics(width=200, height= 100, container=mswin.imgrow3)
			addHandlerChanged(imgmask15, SetLabel15)
			imgmask16 <- ggraphics(width=200, height= 100, container=mswin.imgrow4)
			addHandlerChanged(imgmask16, SetLabel16)
			imgmask17 <- ggraphics(width=200, height= 100, container=mswin.imgrow4)
			addHandlerChanged(imgmask17, SetLabel17)
			imgmask18 <- ggraphics(width=200, height= 100, container=mswin.imgrow4)
			addHandlerChanged(imgmask18, SetLabel18)
			imgmask19 <- ggraphics(width=200, height= 100, container=mswin.imgrow4)
			addHandlerChanged(imgmask19, SetLabel19)
			imgmask20 <- ggraphics(width=200, height= 100, container=mswin.imgrow4)
			addHandlerChanged(imgmask20, SetLabel20)
			imgmask21 <- ggraphics(width=200, height= 100, container=mswin.imgrow5)
			addHandlerChanged(imgmask21, SetLabel21)
			imgmask22 <- ggraphics(width=200, height= 100, container=mswin.imgrow5)
			addHandlerChanged(imgmask22, SetLabel22)
			imgmask23 <- ggraphics(width=200, height= 100, container=mswin.imgrow5)
			addHandlerChanged(imgmask23, SetLabel23)
			imgmask24 <- ggraphics(width=200, height= 100, container=mswin.imgrow5)
			addHandlerChanged(imgmask24, SetLabel24)
			imgmask25 <- ggraphics(width=200, height= 100, container=mswin.imgrow5)
			addHandlerChanged(imgmask25, SetLabel25)
			ref.mswin <- list(imgmask1, imgmask2, imgmask3, imgmask4, imgmask5, imgmask6, imgmask7, imgmask8, imgmask9,
				imgmask10, imgmask11, imgmask12, imgmask13, imgmask14, imgmask15, imgmask16, imgmask17, imgmask18, imgmask19, imgmask20,
				imgmask21, imgmask22, imgmask23, imgmask24, imgmask25)
			
			
			visible(mswin) <- TRUE
			visible(mswin.refimg1) <- TRUE
			par(mar=rep(0,4))
			Image(mimg1[yrange[1]:yrange[2], xrange[1]:xrange[2]], useRaster=T)		
			visible(mswin.refimg2) <- TRUE
			par(mar=rep(0,4))
			Image(mimg2[yrange[1]:yrange[2], xrange[1]:xrange[2]], useRaster=T)
			
	
			# Set the first page to 1
			curpage=1
			# Navigation buttons
			mswin.navgroup <- ggroup(container=mswin.imggroup)
			mswin.prevbutton <- gbutton(text="Previous Page", handler=PagePrev, container=mswin.navgroup)
			enabled(mswin.prevbutton) <- FALSE
			mswin.pagelabel <- glabel(text=paste("Showing 1-", min(25, length(containedmasks)), " of ", length(containedmasks), sep=""), container=mswin.navgroup)
			mswin.nextbutton <- gbutton(text="Next Page", handler=PageNext, container=mswin.navgroup)
			mswin.allnegative <- gbutton(text="Mark unlabeled masks negative", handler=MarkNegative, container=mswin.navgroup)
			mswin.closewindow <- gbutton(text="Close window", handler=CloseLabelWindow, container=mswin.navgroup)
			if(length(containedmasks)<=25){
				enabled(mswin.nextbutton) <- FALSE
			}
			# Plot the first page
			PlotPage(curpage)
		}
	}
	
	
	########## Layout
	
	##The widgets
	win <- gwindow("Calcium Imaging Segmenter", spacing=10)
	maingroup <- ggroup(horizontal=FALSE, cont=win, spacing=5)
	vexpgroup <- ggroup(horizontal=FALSE, cont=maingroup, expand=TRUE, spacing=10)
	controlgroup <- ggroup(horizontal=TRUE, cont=maingroup, expand=TRUE, spacing=10)	
	############################
	##### view experiment group
		
	imgsubgroupA <- ggroup(cont=vexpgroup)
	imgsubgroupB <- ggroup(cont=vexpgroup)

	imgwindowPA <- ggraphics(container=imgsubgroupA)
	imgwindowPB <- ggraphics(container=imgsubgroupB)
	
	imgwindow1 <- ggraphics(container=imgsubgroupA)
	addHandlerChanged(imgwindow1, ImgHandler)
	
	imgwindow2 <- ggraphics(container=imgsubgroupA)
	addHandlerChanged(imgwindow2, ImgHandler)

	imgwindow3 <- ggraphics(container=imgsubgroupB)
	addHandlerChanged(imgwindow3, ImgHandler)

	imgwindow4 <- ggraphics(container=imgsubgroupB)
	addHandlerChanged(imgwindow4, ImgHandler)

	imgwindowPA2 <- ggraphics(container=imgsubgroupA)
	imgwindowPB2 <- ggraphics(container=imgsubgroupB)
	
	visible(imgwindow1) <- TRUE
	par(mar=c(0,0,0,0))
	Image(mimg1, useRaster=T)

	visible(imgwindow2) <- TRUE
	par(mar=c(0,0,0,0))
	Image(mimg2, useRaster=T)

	visible(imgwindow3) <- TRUE
	par(mar=c(0,0,0,0))
	Image(mimg1eq, useRaster=T)

	visible(imgwindow4) <- TRUE
	par(mar=c(0,0,0,0))
	Image(mimg2eq, useRaster=T)

	visible(imgwindowPA) <- TRUE
	par(mar=c(0,0,0,0))
	Image(mimg1, useRaster=T)
	
	visible(imgwindowPA2) <- TRUE
	par(mar=c(0,0,0,0))
	
	Image(mimg2, useRaster=T)
	visible(imgwindowPB) <- TRUE
	par(mar=c(0,0,0,0))
	Image(mimg1eq, useRaster=T)
	
	visible(imgwindowPB2) <- TRUE
	par(mar=c(0,0,0,0))
	Image(mimg2eq, useRaster=T)
	
	
	togglemasks <- gradio(c("No Labels", "Hand Labels", "Final Segmentation"), handler=ToggleMasks, cont=controlgroup)
	
	initiatebutton <- gbutton("Initiate hand labels from classifier", handler=InitiateFromClass, cont=controlgroup)
	redosegbutton <- gbutton("Recompute Segmentation", handler=RedoSeg, cont=controlgroup)

	if(is.null(cf)){
		visible(redosegbutton) <- FALSE
		visible(initiatebutton) <- FALSE
	}
	
	if(all(is.na(selmat[,3]))){
		AddLabeledMasks(F)
		svalue(togglemasks) <- "Hand Labels"
	}else{
		AddLabeledMasks(T)
		svalue(togglemasks) <- "Final Segmentation"
	}
	
	return()
	
}

